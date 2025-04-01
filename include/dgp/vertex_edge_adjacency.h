#ifndef LIBDGP_VERTEX_EDGE_ADJACENCY_H
#define LIBDGP_VERTEX_EDGE_ADJACENCY_H

#include "dgp_inline.h"

#include <Eigen/Dense>

#include <vector>

namespace dgp
{
	template<
		typename DerivedE,
		typename DerivedEMAP,
		typename Index>
	DGP_INLINE void vertex_edge_adjacency(
		const Eigen::PlainObjectBase<DerivedE> &E,
		const Eigen::PlainObjectBase<DerivedEMAP> &EMAP,
		const typename DerivedE::Index nv,
		const typename DerivedE::Index nue,
		std::vector<std::vector<Index> > &VuE);
}

#ifndef DGP_STATIC_LIBRARY
#include   "vertex_edge_adjacency.cpp"
#endif

#endif //LIBDGP_VERTEX_EDGE_ADJACENCY_H
