#ifndef LIBDGP_HALFEDGE_OPPOSITE_H
#define LIBDGP_HALFEDGE_OPPOSITE_H

#include "dgp_inline.h"

#include <Eigen/Dense>

namespace dgp
{
    // input: 
    //  E #F by 3 a mapping from each halfedge to each edge(Unique). 
    //  HeOpp #F by 3 a mapping from each halfedge to its opposite. 
    //  igl::orient_halfedges(F, E, oE);
    // example 1:
    //                       v0         |
    //                     / | \        |
    //                    /  |v3\       |
    //                   /   /\  \      |
    //                  /   /  \ \      |
    //                 / /_ _ _ _\      | 
    //                 v1         v2  
    //
    //  F = [[0, 1, 3], [1, 2, 3], [3, 2, 0]]
    //  E = [[4, 2, 0], [5, 4, 3], [1, 2, 5]]
    //  HeOpp = []
    template <typename DerivedF, typename DerivedHeOpp>
    DGP_INLINE void halfedge_opposite(
        const Eigen::MatrixBase<DerivedF> &F,
        Eigen::PlainObjectBase<DerivedHeOpp> &HeOpp);
}

#ifndef DGP_STATIC_LIBRARY
#include   "halfedge_opposite.cpp"
#endif

#endif //LIBDGP_HALFEDGE_OPPOSITE_H