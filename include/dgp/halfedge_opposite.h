#ifndef LIBDGP_HALFEDGE_OPPOSITE_H
#define LIBDGP_HALFEDGE_OPPOSITE_H

#include "dgp_inline.h"

#include <Eigen/Dense>

namespace dgp
{
    template <typename DerivedE, typename DerivedHeOpp>
    DGP_INLINE void halfedge_opposite(
        const Eigen::MatrixBase<DerivedE> &E,
        Eigen::PlainObjectBase<DerivedHeOpp> &HeOpp);
}

#ifndef DGP_STATIC_LIBRARY
#include   "halfedge_opposite.cpp"
#endif

#endif //LIBDGP_HALFEDGE_OPPOSITE_H