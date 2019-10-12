#ifndef GDB_ARMA_HEADER
#define GDB_ARMA_HEADER
//-------------------------------------------------------------------------
// Supplementary code for debugging Armadillo with gdb, according to Refs.:
// 		https://stackoverflow.com/questions/9499445/is-there-a-way-to-print-an-armadillo-matrix-in-gdb
// 		https://stackoverflow.com/questions/41444181/how-to-print-armadillo-matrix-in-gdb-in-complex-cpp-project?noredirect=1&lq=1
// Usage:
// 		(gdb) call print_matrix<arma::Col<double> >(...)
// The precise command to run in gdb can be found by:
// 		$ nm -C executable | grep print_matrix
//-------------------------------------------------------------------------
#include <armadillo>

template<class Matrix>
void print_matrix(Matrix matrix) {
    matrix.print(std::cout);
}

template void print_matrix<arma::vec>(arma::vec matrix);

#endif // GDB_ARMA_HEADER
