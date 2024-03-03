#ifndef TOOLS_HEADER
#define TOOLS_HEADER
// -------------------------------------------------------
// Copyright (C) 2019 by Dmitry R. Gulevich
// You may use or modify this code but not distribute it.
// ------------------------------------------------------
#include <armadillo>

namespace tdcdft {
namespace tools {

/**
 * Recommend value for OPENBLAS_NUM_THREADS environment variable
 */
void recommend_num_threads(int nthreads);

/**
 * Supplementary code for debugging Armadillo with gdb from the Refs.:
 * https://stackoverflow.com/questions/9499445/is-there-a-way-to-print-an-armadillo-matrix-in-gdb
 * and
 * https://stackoverflow.com/questions/41444181/how-to-print-armadillo-matrix-in-gdb-in-complex-cpp-project?noredirect=1&lq=1
 * To use in gdb debugger call "print_matrix<arma::Col<double> >(...)". The precise command to run can be found by
 * running "nm -C executable | grep print_matrix" in the terminal.
 */
template <class Matrix>
void print_matrix(Matrix matrix) {
    matrix.print(std::cout);
}

template void print_matrix<arma::vec>(arma::vec matrix);

}  // namespace tools
}  // namespace tdcdft

#endif  // TOOLS_HEADER
