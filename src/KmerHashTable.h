#ifndef KALLISTO_KMERHASHTABLE_H
#define KALLISTO_KMERHASHTABLE_H

#include <utility>
#include <string>
#include <iterator>
#include <sys/mman.h>
#include <numa.h>
#include <unistd.h>
#include "common.h"

/*#include <iostream> // debug
	using namespace std;*/

#ifndef MAP_HUGE_1GB
#define MAP_HUGE_1GB (30 << MAP_HUGE_SHIFT)
#endif

#ifndef MAP_HUGE_2MB
#define MAP_HUGE_2MB ((21 << MAP_HUGE_SHIFT))
#endif

// Not thread safe
class KmerHashAllocator {
  private:
  static size_t Roundup(size_t sz, size_t page_sz){
    if (sz & (page_sz-1))
        return (sz & ~(page_sz-1)) + page_sz;
    return sz;
  }
  public:
  static void * Allocate(size_t sz) {
     int flags =  MAP_PRIVATE | MAP_ANONYMOUS;
     int extra_flags[4] = {flags | MAP_HUGETLB | MAP_HUGE_1GB,
                           flags | MAP_HUGETLB | MAP_HUGE_2MB, 
                           flags | MAP_HUGETLB, 
                           flags};
     void * p=MAP_FAILED;

     size_t  default_page_size = sysconf(_SC_PAGESIZE);

     size_t  page_sz[4] = {1<<30,
                       1<<21,
                       default_page_size,
                       default_page_size};
 

     // we hide the size behind the real pointer to use it during deallocation
     sz = sz + sizeof(size_t);
     int cur_flags_idx = 0;
     while (cur_flags_idx < 4) {
         // Round to the nearest page size
         sz = Roundup(sz, page_sz[cur_flags_idx]);
 
         p = mmap(0, sz, PROT_WRITE | PROT_READ, extra_flags[cur_flags_idx], 0, 0);
         if( p != MAP_FAILED) {
             // std::cerr << "[index] Allocated " << sz << " bytes at "<< p <<" attempt #"<<  cur_flags_idx << "\n";
             // remember thr real allocation size
             *(size_t*)p = sz;
             break;
         }
         cur_flags_idx ++;
     }
    
     if (p == MAP_FAILED) {
         // could not allocate the required memory!
         std::cerr << "Error! Failed to mmap() \n";
         exit(-1);
     }

#ifdef NUMA_INTERLEAVE
     // interleave memory on all NUMA nodes, if possible
     if(numa_available() != -1) {
        /* https://linux.die.net/man/2/mbind says:
          Support for huge page policy was added with 2.6.16. 
          For interleave policy to be effective on huge page mappings the policied memory needs to be tens of megabytes or larger.
          If the allocated size is tiny, i.e., < 128 MB, do not bother interleaving. Instead take advantage of the caching on a processor.
        */
        if( sz >= (1L<<27)) {
            // We call numa_get_mems_allowed() each time to get the possibly changing list
            struct bitmask * mem_nodes = numa_get_mems_allowed();
            numa_interleave_memory(p, sz, mem_nodes);
        }
     }
#endif
     p = (void *)((size_t)p + sizeof(size_t));
     return p;
 }

  static void Deallocate(void * p){
      // The real pointer is size_t behind p.
      p = (void*)((size_t)p - sizeof(size_t));
      // the real size is stored one behind
      size_t sz = * (size_t*) p;
      //printf("\nDeallocated  %lu at %p \n", sz, p);
      if (munmap(p, sz) != 0 ) {
          std::cerr << "Error! Failed to munmap() \n";
          exit(-1);
      }
  }
};



template<typename T, typename Hash = KmerHash>
struct KmerHashTable {
  using value_type = std::pair<Kmer, T>;
  using key_type = Kmer;
  using mapped_type = T;

  Hash hasher;
  value_type *table;
  size_t size_, pop;
  value_type empty;
  value_type deleted;
  double load_factor;
  size_t threshold_to_resize;


// ---- iterator ----

  template<bool is_const_iterator = true>
  class iterator_ : public std::iterator<std::forward_iterator_tag, value_type> {
   public:

    typedef typename std::conditional<is_const_iterator, const KmerHashTable *, KmerHashTable *>::type DataStructurePointerType;
    typedef typename std::conditional<is_const_iterator, const value_type&, value_type&>::type ValueReferenceType;
    typedef typename std::conditional<is_const_iterator, const value_type *, value_type *>::type ValuePointerType;


    DataStructurePointerType ht;
    size_t h;

    iterator_() : ht(nullptr), h(0) {}
    iterator_(DataStructurePointerType ht_) : ht(ht_), h(ht_->size_) {}
    iterator_(DataStructurePointerType ht_, size_t h_) :  ht(ht_), h(h_) {}
    iterator_(const iterator_<false>& o) : ht(o.ht), h(o.h) {}
    iterator_& operator=(const iterator_& o) {ht=o.ht; h=o.h; return *this;}

    ValueReferenceType operator*() const {return ht->table[h];}
    ValuePointerType operator->() const {return &(ht->table[h]);}

    void find_first() {
      h = 0;
      if (ht->table != nullptr && ht->size_>0) {
        Kmer& km = ht->table[h].first;
        if (km == ht->empty.first || km == ht->deleted.first) {
          operator++();
        }
      }
    }

    iterator_ operator++(int) {
      const iterator_ old(*this);
      ++(*this);
      return old;
    }

    iterator_& operator++() {
      if (h == ht->size_) {
        return *this;
      }
      ++h;
      for (; h < ht->size_; ++h) {
        Kmer& km = ht->table[h].first;
        if (km != ht->empty.first && km != ht->deleted.first) {
          break;
        }
      }
      return *this;
    }
    bool operator==(const iterator_ &o) const {return (ht->table == o.ht->table) && (h == o.h);}
    bool operator!=(const iterator_ &o) const {return !(this->operator==(o));}
    friend class iterator_<true>;
  };

  typedef iterator_<true> const_iterator;
  typedef iterator_<false> iterator;


  // --- hash table


  KmerHashTable(const double ht_load_factor = DEFAULT_HT_LOAD_FACTOR, const Hash& h = Hash()) : load_factor(ht_load_factor), hasher(h), table(nullptr), size_(0), pop(0) {
    empty.first.set_empty();
    deleted.first.set_deleted();
    init_table(1024);
  }

  KmerHashTable(size_t sz, const double ht_load_factor = DEFAULT_HT_LOAD_FACTOR, const Hash& h = Hash()) : hasher(h), table(nullptr), size_(0), load_factor(ht_load_factor), pop(0)  {
    empty.first.set_empty();
    deleted.first.set_deleted();
    // init a table of size that is slightly larger than sz*1.0/ht_load_factor
    init_table((size_t) (sz*1.0/ht_load_factor + (1L<<20)));
  }

  ~KmerHashTable() {
    clear_table();
  }

  void clear_table() {
    if (table != nullptr) {
#ifdef USE_CUSTOM_HASH_ALLOCATOR
      //  delete[] table;
      // We cannot call delete because it was allocated with a c++ placement allocator
      table->~value_type(); // explicit destrictor
      KmerHashAllocator::Deallocate(table); // deallocate
#else
     delete[] table;
#endif
      table = nullptr;
    }
    size_ = 0;
    pop  = 0;
  }

  size_t size() const {
    return pop;
  }

  void clear() {
    std::fill(table, table+size_, empty);
    pop = 0;
  }

  void init_table(size_t sz) {
    clear_table();
    size_ = rndup(sz);
    threshold_to_resize = size_ * load_factor;
    //cerr << "init table of size " << size_ << endl;
    // A placement allocator with MMAP and huge TLB pages
#ifdef USE_CUSTOM_HASH_ALLOCATOR
    void * ptr = KmerHashAllocator::Allocate(size_ * sizeof(value_type));
    table = new (ptr) value_type[size_];
    // No need to zero init 
#else
    table = new value_type[size_];
    std::fill(table, table+size_, empty);
#endif
  }

  iterator find(const Kmer& key) {
    size_t h = hasher(key) & (size_-1);

    for (;; h =  (h+1!=size_ ? h+1 : 0)) {
      if (table[h].first == empty.first) {
        // empty slot, not in table
        return iterator(this);
      } else if (table[h].first == key) {
        // same key, found
        return iterator(this, h);
      } // if it is deleted, we still have to continue
    }
  }

  const_iterator find(const Kmer& key) const {

    size_t h = hasher(key) & (size_-1);

    for (;; h =  (h+1!=size_ ? h+1 : 0)) {
      if (table[h].first == empty.first) {
        // empty slot, not in table
        return const_iterator(this);
      } else if (table[h].first == key) {
        // same key, found
        return const_iterator(this, h);
      }
    }
  }


  iterator erase(const_iterator pos) {
    if (pos == this->end()) {
      return this->end();
    }
    size_t h = pos.h;
    table[h] = deleted;
    --pop;
    return ++iterator(this, h); // return pointer to next element
  }

  size_t erase(const Kmer& km) {
    const_iterator pos = find(km);
    size_t oldpop = pop;
    if (pos != this->end()) {
      erase(pos);
    }
    return oldpop-pop;
  }

  std::pair<iterator,bool> insert(const value_type& val) {
    //cerr << "inserting " << val.first.toString() << " = " << val.second << endl;
    if (pop >= threshold_to_resize) {
      //cerr << "-- triggered resize--" << endl;
      reserve(2*size_);
    }

    size_t h = hasher(val.first) & (size_-1);
    //cerr << " hash value = " << h << endl;
    for (;; h = (h+1!=size_ ? h+1 : 0)) {
      //cerr << "  lookup at " << h << endl;
      if (table[h].first == empty.first || table[h].first == deleted.first) {
        //cerr << "   found empty slot" << endl;
        // empty slot, insert here
        table[h] = val;
        ++pop; // new table
        return {iterator(this, h), true};
      } else if (table[h].first == val.first) {
        // same key, update value
        //cerr << "   found key already here " << table[h].first.toString() << " = " << table[h].second <<  endl;
        return {iterator(this, h), false};
      }
    }

  }

  void reserve(size_t sz) {
    if (sz <= size_) {
      return;
    }

    value_type *old_table = table;
    size_t old_size_ = size_;


    size_ = rndup(sz);
    threshold_to_resize = size_ * load_factor;
    pop = 0;
#ifdef USE_CUSTOM_HASH_ALLOCATOR
    // A placement allocator with MMAP and huge TLB pages
    void * ptr = KmerHashAllocator::Allocate(size_ * sizeof(value_type));
    table = new (ptr) value_type[size_];
    // no need to init the table since MMAP gives zerored out pages
#else
    table = new value_type[size_];
    std::fill(table, table+size_, empty);
#endif

    for (size_t i = 0; i < old_size_; i++) {
      if (old_table[i].first != empty.first && old_table[i].first != deleted.first) {
        insert(old_table[i]);
      }
    }

#ifdef USE_CUSTOM_HASH_ALLOCATOR
    // We cannot call delete because it was allocated with a c++ placement allocator
    old_table->~value_type(); // explicit destrictor
    KmerHashAllocator::Deallocate(old_table); // deallocate
#else
    delete[] old_table;
#endif
    old_table = nullptr;
  }

  size_t rndup(size_t v) {
    v--;
    v |= v >> 1;
    v |= v >> 2;
    v |= v >> 4;
    v |= v >> 8;
    v |= v >> 16;
    v |= v >> 32;
    v++;
    return v;
  }

  iterator begin() {
    iterator it(this);
    it.find_first();
    return it;
  }

  const_iterator begin() const {
    const_iterator it(this);
    it.find_first();
    return it;
  }

  iterator end() {
    return iterator(this);
  }

  const_iterator end() const {
    return const_iterator(this);
  }





};

#endif // KALLISTO_KMERHASHTABLE_H
