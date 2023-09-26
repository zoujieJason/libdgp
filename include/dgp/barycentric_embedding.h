#ifndef LIBDGP_BARYCENTRIC_EMBEDDING_H
#define LIBDGP_BARYCENTRIC_EMBEDDING_H

#include "dgp_inline.h"

#include <Eigen/Dense>

#include <vector>

namespace dgp
{
    DGP_INLINE bool barycentric_embedding(
        const Eigen::MatrixXd &V,
        const Eigen::MatrixXi &F,
        const std::vector<size_t> &mesh_boundary,
        const Eigen::MatrixXd &embedding_boundary,
        Eigen::MatrixXd &UV);
}

#ifndef DGP_STATIC_LIBRARY
#include   "barycentric_embedding.cpp"
#endif

#endif //LIBDGP_BARYCENTRIC_EMBEDDING_H