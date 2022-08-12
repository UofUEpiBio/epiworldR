#include "cpp11.hpp"
using namespace cpp11;
namespace writable = cpp11::writable;

//' This thing sums vectors
//' @param x A numeric vector
//' @export
[[cpp11::register]]
double sum_cpp(doubles x) {
  int n = x.size();
  double sum = 0;
  for(int i = 0; i < n; ++i) {
    sum += x[i];
  }
  return sum;
}

// cpp11::cpp_source("~/Desktop/Grad Research/sum_function.cpp")