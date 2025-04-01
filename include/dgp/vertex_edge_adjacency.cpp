#include "vertex_edge_adjacency.h"

namespace dgp
{
	template<
		typename DerivedE,
		typename DerivedEMAP,
		typename Index>
	void vertex_edge_adjacency(
		const Eigen::PlainObjectBase<DerivedE>& E,
		const Eigen::PlainObjectBase<DerivedEMAP>& EMAP,
		const typename DerivedE::Index nv,
		const typename DerivedE::Index nue,
		std::vector<std::vector<Index> >& VuE)
	{
		VuE.resize(nv);
		std::vector<std::vector<Index> > VE(nv);
		for (Index eid = 0; eid < E.rows(); ++eid)
		{
			VE[E(eid, 0)].push_back(eid);
			VE[E(eid, 1)].push_back(eid);
		}
		for (Index vid = 0; vid < nv; ++vid)
		{
			std::vector<bool> uE_vis(nue, false);
			for (const auto &eid : VE[vid])
			{
				const auto ueid = EMAP(eid);
				if (uE_vis[ueid])
				{
					continue;
				}
				uE_vis[ueid] = true;
				VuE[vid].push_back(ueid);
			}
		}
	}
}