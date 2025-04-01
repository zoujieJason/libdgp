#ifndef LIBDGP_FILLING_HOLES_H
#define LIBDGP_FILLING_HOLES_H

#include "dgp_inline.h"

#include <Eigen/Dense>

#include <vector>

namespace dgp
{
	
	template <
		typename DerivedV,
		typename DerivedF,
		typename DerivedVr,
		typename DerivedFr>
	DGP_INLINE bool filling_holes(
		const Eigen::PlainObjectBase<DerivedV> &V,
		const Eigen::PlainObjectBase<DerivedF> &F,
		Eigen::PlainObjectBase<DerivedVr> &Vr,
		Eigen::PlainObjectBase<DerivedFr> &Fr);

	template <
		typename DerivedV,
		typename DerivedF,
		typename DerivedFr>
	DGP_INLINE bool filling_holes(
		const Eigen::PlainObjectBase<DerivedV> &V,
		const Eigen::PlainObjectBase<DerivedF> &F,
		std::vector<int> &polyline_of_hole,
		Eigen::PlainObjectBase<DerivedFr> &Fr);
}

#ifndef DGP_STATIC_LIBRARY
#include   "filling_holes.cpp"
#endif

#endif //LIBDGP_FILLING_HOLES_H