#ifndef KALLISTO_COMMON_H
#define KALLISTO_COMMON_H

#define KALLISTO_VERSION "0.43.1"

#include <string>
#include <vector>

#define DEFAULT_HT_LOAD_FACTOR (0.8)

struct ProgramOptions {
  bool verbose;
  int threads;
  std::string index;
  int k;
  int iterations;
  std::string output;
  std::string output_filename_prefix; 
  int skip;
  size_t seed;
  double fld;
  double sd;
  double ht_load_factor;
  int min_range;
  int bootstrap;
  std::vector<std::string> transfasta;
  bool batch_mode;
  std::string batch_file_name;
  std::vector<std::vector<std::string>> batch_files;
  std::vector<std::string> batch_ids;
  std::vector<std::string> files;
  std::vector<std::string> umi_files;
  bool plaintext;
  bool write_index;
  bool single_end;
  bool strand_specific;
  bool peek; // only used for H5Dump
  bool bias;
  bool pseudobam;
  bool make_unique;
  bool estimated_counts;
  bool fusion;
  enum class StrandType {None, FR, RF};
  StrandType strand;
  bool umi;
  std::string gfa; // used for inspect

ProgramOptions() :
  verbose(false),
  threads(1),
  k(31),
  iterations(500),
  skip(1),
  seed(42),
  fld(0.0),
  sd(0.0),
  min_range(1),
  bootstrap(0),
  batch_mode(false),
  plaintext(false),
  write_index(false),
  single_end(false),
  strand_specific(false),
  peek(false),
  bias(false),
  pseudobam(false),
  make_unique(false),
  fusion(false),
  strand(StrandType::None),
  umi(false),
  estimated_counts(false),
  output_filename_prefix(""),
  ht_load_factor(DEFAULT_HT_LOAD_FACTOR)
  {}
};

std::string pretty_num(size_t num);
std::string pretty_num(int num);

#endif // KALLISTO_COMMON_H
