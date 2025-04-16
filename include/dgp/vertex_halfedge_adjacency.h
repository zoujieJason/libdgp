#ifndef LIBDGP_VERTEX_HALFEDGE_ADJACENCY_H
#define LIBDGP_VERTEX_HALFEDGE_ADJACENCY_H

#include "dgp_inline.h"

#include <Eigen/Dense>

#include <tuple>

namespace dgp
{
    // usage: 
    //      for(auto vhe: VHe[v])
    //      {
    //          int f = std::get<0>(vhe);
    //          int k = std::get<1>(vhe);   // k -> e_ij; outgoing he;
    //          bool outgoing = std::get<3>(vhe);
    //          int i = (k + 1) % 3; // i -> e_jk; f(i,0) = v;
    //      }
    //      const auto fi2e = [&nf](int f, int i) ->int { return nf * i + f; }; 
    //      const auto e2f = [&nf](int e)->int { return e % nf; };
    //      const auto e2fi = [&nf](int e)->int { return e / nf; }; // vertex F(f, i) opposite to edge e
    // example 1: 
    //  F = [ [0, 1, 2] ]         
    //  HeOpp = [-1, -1, -1]  // see halfedge_opposite
    // 
    //  In this case:
    //  VHe[0] = [[0, 2]]
    //  VHe[1] = [[0, 0]]
    //  VHe[2] = [[0, 1]]

    DGP_INLINE void decode_halfedge(
        const std::tuple<int, int, bool> &he,
        int &f, int &i, int &j, int &k);

    DGP_INLINE int vertex_traversal_find_ingoing(
        const Eigen::MatrixXi &F,
        const Eigen::MatrixXi &HeOpp,
        int v, 
        int outgoing_he, 
        std::vector<std::tuple<int, int, bool>> &helist,
        int &ingoing_he);

    DGP_INLINE int vertex_traversal_find_outgoing(
        const Eigen::MatrixXi &F,
        const Eigen::MatrixXi &HeOpp,
        int v, 
        int ingoing_he, 
        std::vector<std::tuple<int, int, bool>> &helist,
        int &outgoing_he);

    DGP_INLINE bool vertex_halfedge_adjacency(
        const std::vector<std::vector<int> > &VF, 
        const std::vector<std::vector<int> > &VFi,
        const Eigen::MatrixXi &F,
        const Eigen::MatrixXi &HeOpp,
        std::vector<std::pair<bool, std::vector<std::tuple<int, int, bool>>>> &VHe);
}

#ifndef DGP_STATIC_LIBRARY
#include   "vertex_halfedge_adjacency.cpp"
#endif

#endif //LIBDGP_VERTEX_HALFEDGE_ADJACENCY_H