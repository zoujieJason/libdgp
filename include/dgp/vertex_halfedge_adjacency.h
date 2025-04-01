#ifndef LIBDGP_VERTEX_HALFEDGE_ADJACENCY_H
#define LIBDGP_VERTEX_HALFEDGE_ADJACENCY_H

#include "dgp_inline.h"

#include <Eigen/Dense>

namespace dgp
{
    DGP_INLINE bool vertex_halfedge_adjacency(
        const std::vector<std::vector<int> > &VF, 
        const std::vector<std::vector<int> > &VFi,
        const Eigen::MatrixXi &F,
        const Eigen::MatrixXi &HeOpp,
        const std::vector<bool> &is_border_vertex,
        std::vector<std::vector<std::pair<int, int>>> &VHe);
}

#ifndef DGP_STATIC_LIBRARY
#include   "vertex_halfedge_adjacency.cpp"
#endif

#endif //LIBDGP_VERTEX_HALFEDGE_ADJACENCY_H