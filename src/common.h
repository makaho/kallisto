#ifndef KALLISTO_COMMON_H
#define KALLISTO_COMMON_H

#define KALLISTO_VERSION "0.42"

#include <string>
#include <vector>

struct ProgramOptions {
  bool verbose;
  int threads;
  std::string index;
  int k;
  int iterations;
  std::string output;
  int skip;
  size_t seed;
  double fld;
  int min_range;
  int bootstrap;
  std::string transfasta;
  std::vector<std::string> files;
  bool plaintext;

ProgramOptions() :
  verbose(false),
  seed(42),
  threads(1),
  k(31),
  iterations(500),
  skip(1),
  min_range(1),
  fld(0.0),
  bootstrap(0),
  plaintext(false)
  {}
};

#endif // KALLISTO_COMMON_H
