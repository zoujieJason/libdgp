#include "vertex_halfedge_adjacency.h"

#include <igl/is_border_vertex.h>

namespace dgp
{
    DGP_INLINE bool vertex_halfedge_adjacency(
        const std::vector<std::vector<int> > &VF, 
        const std::vector<std::vector<int> > &VFi,
        const Eigen::MatrixXi &F,
        const Eigen::MatrixXi &HeOpp,
        const std::vector<bool> &is_border_vertex,
        std::vector<std::vector<std::pair<int, int>>> &VHe)
    {
        // he: a first outgoing he of vertex v
        const auto ccw = [&](
            int v, 
            int he,
            std::vector<std::pair<int, int>> &helist,
            int &outgoing_he)->bool 
        {
            helist.clear();
            outgoing_he = he;
            while(outgoing_he != -1)
            {
                int f = outgoing_he % F.rows();
                int k = outgoing_he / F.rows();
                helist.push_back(std::make_pair(f, k));
                int i = (k + 1) % 3;
                if(F(f, i) != v) return false;
                int j = (k + 2) % 3;
                outgoing_he = HeOpp(f, j);
            } 
            return true; 
        };

        for(size_t v = 0; v < VF.size(); ++v)
        {
            if(!VF[v].empty()) return false;
            
            int f0 = VF[v].front(); 
            int i0 = VFi[v].front(); // i -> e(j, k);
            int j0 = (i0 + 1) % 3;   // j -> e(k, i);
            int k0 = (i0 + 2) % 3;   // k -> e(i, j);

            std::vector<std::pair<int, int>> helist;
            int outgoing_he; 
            if(!ccw(v, f0 + k0 * F.rows(), helist, outgoing_he)) return false; 

            if(outgoing_he == -1) 
            {
                if(!is_border_vertex[v]) return false;
                helist.clear(); 
                outgoing_he = f0 + k0 * F.rows();
                int border_he = outgoing_he;
                outgoing_he = HeOpp(f0, k0); 
                while(outgoing_he != -1)
                {
                    int f = outgoing_he % F.rows();
                    int j = outgoing_he / F.rows();
                    int k = (j + 1) % 3;
                    int i = (j + 2) % 3;
                    if(F(f, i) != v) return false;
                    border_he = f + k * F.rows();
                    outgoing_he = HeOpp(f, k);
                }
                if(!ccw(v, border_he, helist, outgoing_he)) return false; 
                if(outgoing_he != -1) return false;
            }
            VHe.push_back(helist);
        }
        return true; 
    }
}

#ifdef DGP_STATIC_LIBRARY

#endif
