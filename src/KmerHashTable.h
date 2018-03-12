#ifndef KALLISTO_KMERHASHTABLE_H
#define KALLISTO_KMERHASHTABLE_H

#include <utility>
#include <string>
#include <iterator>
#include <sys/mman.h>
#include <numa.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdint.h>
#include <limits.h>
#include "common.h"

#define PROFILE_CLASH

#ifdef PROFILE_CLASH
	extern uint64_t gLoop;	
	extern uint64_t gMaxLoop;
	extern uint64_t gInvoke;
#endif
 

/*#include <iostream> // debug
	using namespace std;*/
using namespace std;
#ifndef MAP_HUGE_1GB
#define MAP_HUGE_1GB (30 << MAP_HUGE_SHIFT)
#endif

#ifndef MAP_HUGE_2MB
#define MAP_HUGE_2MB ((21 << MAP_HUGE_SHIFT))
#endif

union MapMetadata {
    struct {
        size_t allocSz;
        size_t serPop;
        size_t serSize;
        int fd;
    };
    char dummy[128];
};

#define MAPDATA_TO_TABLE(_p) ((void *)((size_t)_p + sizeof(MapMetadata)))
#define TABLE_TO_MAPDATA(_p) ((MapMetadata*)((size_t)_p - sizeof(MapMetadata)))

// Not thread safe
class KmerHashAllocator {
  private:
  size_t Roundup(size_t sz, size_t page_sz){
    if (sz & (page_sz-1))
        return (sz & ~(page_sz-1)) + page_sz;
    return sz;
  }
    
    void * AllocWithHugePage(int flags, const size_t sz, const int protection, const int fd, const bool updateMetadata=true){
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
        int cur_flags_idx = 0;
        while (cur_flags_idx < 4) {
            // Round to the nearest page size
            size_t realSz = Roundup(sz, page_sz[cur_flags_idx]);
            
            p = mmap(0, realSz, protection, extra_flags[cur_flags_idx], fd, 0);
            if( p != MAP_FAILED) {
                std::cerr << "[index] Allocated " << realSz << " bytes at "<< p <<" attempt #"<<  cur_flags_idx << "\n";
                // remember thr real allocation size
                if (updateMetadata == true) {
                    MapMetadata* metaData = (MapMetadata*)p;
                    metaData->allocSz = realSz;
                    metaData->fd = fd;
                }
                break;
            }
            cur_flags_idx ++;
        }
        return p;
    }
  public:
  void * Allocate(size_t sz) {
     int flags =  MAP_PRIVATE | MAP_ANONYMOUS;
     int protection = PROT_WRITE | PROT_READ;
     // we hide the size behind the real pointer to use it during deallocation
     sz = sz + sizeof(MapMetadata);
     void * p = AllocWithHugePage(flags, sz, protection, 0 /*fd*/);
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
	size_t realsz = ((MapMetadata*)(p))->allocSz;
        if( sz >= (1L<<27)) {
            // We call numa_get_mems_allowed() each time to get the possibly changing list
            struct bitmask * mem_nodes = numa_get_mems_allowed();
            numa_interleave_memory(p, realsz, mem_nodes);
        }
     }
#endif
     return MAPDATA_TO_TABLE(p);
 }

  void Deallocate(void * p){
      // Metadata is behind p.
      MapMetadata * metaData = TABLE_TO_MAPDATA(p);
      if (munmap(metaData, metaData->allocSz) != 0 ) {
          std::cerr << "Error! Failed to munmap() \n";
          exit(-1);
      }
  }

  void WarnHugeTLB(const string &file) {
	  char path[2*PATH_MAX];
	  char canonicalPath[2*PATH_MAX];
	  char * saveptr;
	  char * mountpt;
	  const char * delim(" ");
	  /* Open the command for reading. */
	  FILE * fp = popen("mount -t hugetlbfs", "r");
	  if (fp == NULL) 
		  goto ErrRet;
	  // output: hugetlbfs on /dev/hugepages type hugetlbfs (rw,relatime)
	  if ( NULL == fgets(path, 10*PATH_MAX-1, fp)) {
		  goto ErrRet;
	  }
	  pclose(fp);
	  fp = NULL;
	  if (0 == strtok_r(path, delim, &saveptr))
		  goto ErrRet;
	  if ( 0 == strtok_r(NULL, delim, &saveptr))
		  goto ErrRet;
	  if ( (mountpt=strtok_r(NULL, delim, &saveptr)) == 0 )
		  goto ErrRet;
	  if (NULL == realpath(mountpt, canonicalPath))
		  goto ErrRet;
	  if (NULL == strstr(file.c_str(), canonicalPath))
		  std::cerr << "[Warn]" << file << " does not seem to be on hugetlbfs. It can severly affect performance" << endl;
	  return;
ErrRet:
	  if (fp)
		  pclose(fp);
	  std::cerr << "[Warn] Could not verify if " << file << " is on hugetlbfs" << endl;
  }

  bool Serialize(void * p, string file, size_t __size, size_t __pop){
	  WarnHugeTLB(file);
	  // Interleave the allocation
#ifdef NUMA_INTERLEAVE
	  // interleave memory on all NUMA nodes, if possible
	  if(numa_available() != -1) {
		  struct bitmask * mem_nodes = numa_get_mems_allowed();
		  numa_set_interleave_mask(mem_nodes);
	  }
#endif
	  int fd = open(file.c_str(), O_CREAT|O_TRUNC|O_RDWR, S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP);
	  if (fd == -1) {
		  perror("open");
		  return false;
	  }
	  // Metadata is behind p.
	  MapMetadata * metaData = TABLE_TO_MAPDATA(p);
	  metaData->serSize = __size;
	  metaData->serPop = __pop;
	  int ret = ftruncate(fd, metaData->allocSz);
	  if (ret == -1) {
		  perror("ftruncate");
		  return false;
	  }

	  int flags =  MAP_SHARED;
	  int protection = PROT_WRITE | PROT_READ;
	  void * fileBackedData = AllocWithHugePage(flags, metaData->allocSz, protection, fd);
	  if (fileBackedData == MAP_FAILED) {
		  // could not allocate the required memory!
		  std::cerr << "Error! Failed to mmap() \n";
		  return false;
	  }
	  MapMetadata * fileBackedMetaData = (MapMetadata*)(fileBackedData);
	  size_t serializedAllocSz = fileBackedMetaData->allocSz;	
	  // sanity check that fileBackedMetaData->allocSz >= metaData->allocSz
	  assert (serializedAllocSz >= metaData->allocSz);
	  memcpy(fileBackedData, metaData, metaData->allocSz);
	  // memcpy overwrites fileBackedMetaData->allocSz, hence reset it
	  fileBackedMetaData->allocSz = serializedAllocSz;
	  if (close(fd) == -1) {
		  perror("close");
		  return false;
	  }
	  // This guarantees that all data is written back.
	  if (munmap(fileBackedData, metaData->allocSz) != 0 ) {
		  std::cerr << "Error! Failed to munmap() \n";
		  exit(-1);
	  }
	  return true;
  }

  void * Deserialize(string file, size_t & ref_size, size_t & ref_pop){
	  WarnHugeTLB(file);
	  int fd = open(file.c_str(), O_RDONLY);
	  if (fd == -1) {
		  perror("open");
		  return 0;
	  }

	  // get file size
	  struct stat stFileInfo;
	  auto intStat = stat(file.c_str(), &stFileInfo);
	  if (intStat != 0) {
		  throw std::runtime_error ("Cannot stat" + file);
	  }
	  size_t sz = stFileInfo.st_size;

	  int flags =  MAP_SHARED;
	  int protection = PROT_READ;
	  void * fileBackedData = AllocWithHugePage(flags, sz, protection, fd, false /*updateMetadata*/);
	  if (fileBackedData == MAP_FAILED) {
		  // could not allocate the required memory!
		  throw std::runtime_error ("Cannot memory map" + file);
	  }
	  // Close immediately.
	  // ref: http://pubs.opengroup.org/onlinepubs/7908799/xsh/mmap.html
	  // The mmap() function adds an extra reference to the file associated
	  // with the file descriptor fildes which is not removed by a subsequent close()
	  // on that file descriptor. This reference is removed when there are no more mappings to the file.
	  if (close(fd) == -1) {
		  perror("close");
		  // silently continue
	  }

	  // sanity check
	  MapMetadata * metaData = (MapMetadata*)(fileBackedData);
	  if(metaData->allocSz != sz) {
		  throw std::runtime_error ("Something went wrong between serialization and deserialization of " + file);
	  }

	  ref_size = metaData->serSize;
	  ref_pop =  metaData->serPop;
	  return MAPDATA_TO_TABLE(fileBackedData);
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
  KmerHashAllocator hashAllocator;


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

    KmerHashTable(string file, const Hash& h = Hash()) : hasher(h), table(nullptr), size_(0), pop(0)  {
        empty.first.set_empty();
        deleted.first.set_deleted();
        hashAllocator.Deserialize(file, size_, pop);
    }

  ~KmerHashTable() {
    clear_table();
  }

  void clear_table() {
    if (table != nullptr) {
#ifdef USE_CUSTOM_HASH_ALLOCATOR
      //  delete[] table;
      // We cannot call delete because it was allocated with a c++ placement allocator
      table->~value_type(); // explicit destrctor
      hashAllocator.Deallocate(table); // deallocate
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
    void * ptr = hashAllocator.Allocate(size_ * sizeof(value_type));
    table = new (ptr) value_type[size_];
    // No need to zero init 
#else
    table = new value_type[size_];
    std::fill(table, table+size_, empty);
#endif
  }

  iterator find(const Kmer& key) {
    size_t h = hasher(key) & (size_-1);
#ifdef PROFILE_CLASH
    	int loopTrip = 0;
	gInvoke++;
#endif


    for (;; h =  (h+1!=size_ ? h+1 : 0)) {
#ifdef PROFILE_CLASH
	loopTrip++;	
#endif
      if (table[h].first == empty.first) {
#ifdef PROFILE_CLASH
	gLoop += loopTrip;	
	if (loopTrip > gMaxLoop)
		gMaxLoop = loopTrip;
#endif
        // empty slot, not in table
        return iterator(this);
      } else if (table[h].first == key) {
#ifdef PROFILE_CLASH
	gLoop += loopTrip;	
	if (loopTrip > gMaxLoop)
		gMaxLoop = loopTrip;
#endif
        // same key, found
        return iterator(this, h);
      } // if it is deleted, we still have to continue
    }
  }

  const_iterator find(const Kmer& key) const {

    size_t h = hasher(key) & (size_-1);
#ifdef PROFILE_CLASH
    	int loopTrip = 0;
	gInvoke++;
#endif


    for (;; h =  (h+1!=size_ ? h+1 : 0)) {
#ifdef PROFILE_CLASH
	loopTrip++;	
#endif
 
      if (table[h].first == empty.first) {
#ifdef PROFILE_CLASH
	gLoop += loopTrip;	
	if (loopTrip > gMaxLoop)
		gMaxLoop = loopTrip;
#endif
 
        // empty slot, not in table
        return const_iterator(this);
      } else if (table[h].first == key) {
#ifdef PROFILE_CLASH
	gLoop += loopTrip;	
	if (loopTrip > gMaxLoop)
		gMaxLoop = loopTrip;
#endif
 
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
    void * ptr = hashAllocator.Allocate(size_ * sizeof(value_type));
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
    hashAllocator.Deallocate(old_table); // deallocate
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

  // -- serialization
  bool write_to_file(string file) {
      return hashAllocator.Serialize(table, file, size_, pop);
  }
    
  bool load_from_file(string file) {
    clear_table();
    table = (value_type*) hashAllocator.Deserialize(file, size_, pop);
    return true;
  }

};

#endif // KALLISTO_KMERHASHTABLE_H
