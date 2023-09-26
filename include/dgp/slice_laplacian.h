#ifndef LIBDGP_SLICE_H
#define LIBDGP_SLICE_H

#include "dgp_inline.h"

#include <Eigen/Dense>

#include <vector>

namespace dgp
{
    DGP_INLINE bool slice_laplacian(
        const Eigen::SparseMatrix<double> &L,
        const std::vector<size_t> &boundary,
        Eigen::SparseMatrix<double> &LA,
        Eigen::SparseMatrix<double> &LB,
        std::vector<bool> &is_boundary, 
        std::vector<size_t> &vertices_map);
}

#ifndef DGP_STATIC_LIBRARY
#include   "slice_laplacian.cpp"
#endif

#endif //LIBDGP_SLICE_H