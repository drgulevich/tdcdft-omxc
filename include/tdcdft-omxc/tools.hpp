#ifndef TOOLS_HEADER
#define TOOLS_HEADER

#include <armadillo>

/**-------------------------------------------------------------------------
* Supplementary code for debugging Armadillo with gdb, according to Refs.:
* 		https://stackoverflow.com/questions/9499445/is-there-a-way-to-print-an-armadillo-matrix-in-gdb
* 		https://stackoverflow.com/questions/41444181/how-to-print-armadillo-matrix-in-gdb-in-complex-cpp-project?noredirect=1&lq=1
* Usage:
* 		(gdb) call print_matrix<arma::Col<double> >(...)
* The precise command to run in gdb can be found by:
* 		$ nm -C executable | grep print_matrix
*-------------------------------------------------------------------------*/
template<class Matrix>
void print_matrix(Matrix matrix) {
    matrix.print(std::cout);
}

template void print_matrix<arma::vec>(arma::vec matrix);

/** 
* Recommend value for OPENBLAS_NUM_THREADS environment variable
*/
void recommend_num_threads(int nthreads);

#endif // TOOLS_HEADER
