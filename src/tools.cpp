// -------------------------------------------------------
// Copyright (C) 2019 by Dmitry R. Gulevich
// You may use or modify this code but not distribute it.
// -------------------------------------------------------
#include <stdlib.h>  // getenv

#include <iostream>

namespace tdcdft {
namespace tools {

// Recommend value for OPENBLAS_NUM_THREADS environment variable
void recommend_num_threads(int nthreads) {
    char *envvar = getenv("OPENBLAS_NUM_THREADS");
    if (envvar == NULL) std::cout << "# OPENBLAS_NUM_THREADS is not set" << std::endl;
    std::cout << "# export OPENBLAS_NUM_THREADS=" << nthreads << " is recommended for this problem" << std::endl;
}

}  // namespace tools
}  // namespace tdcdft
