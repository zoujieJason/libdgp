#include "vertex_halfedge_adjacency.h"

namespace dgp
{
    DGP_INLINE void decode_halfedge(const std::tuple<int, int, bool> &he, int &f, int &i, int &j, int &k)
    {
        f = std::get<0>(he);
        if(std::get<2>(he))
        {
            k = std::get<1>(he);
            i = (k + 1) % 3;
            j = (i + 1) % 3; 
        }
        else 
        {
            j = std::get<1>(he);
            k = (j + 1) % 3;
            i = (k + 1) % 3;
        }
    }

    DGP_INLINE int vertex_traversal_find_ingoing(
        const Eigen::MatrixXi &F,
        const Eigen::MatrixXi &HeOpp,
        int v, 
        int outgoing_he, 
        std::vector<std::tuple<int, int, bool>> &helist,
        int &ingoing_he)
    {
        const auto nf = F.rows();
        helist.clear();

        int f, k, i, j;
        int oHe = outgoing_he;
        do
        {
            f = oHe % nf;
            k = oHe / nf;
            i = (k + 1) % 3;
            j = (i + 1) % 3;
            if(F(f, i) != v)
            {
                return -1;
            }

            helist.emplace_back(std::make_tuple(f, k, true));
            ingoing_he = f + j * nf;
            oHe = HeOpp(f, j);
            if(oHe == outgoing_he)
            {
                return 1;
            }
        } while (oHe != -1);
        helist.emplace_back(std::make_tuple( ingoing_he % nf, ingoing_he / nf, false));
        return 0; 
    }

    DGP_INLINE int vertex_traversal_find_outgoing(
        const Eigen::MatrixXi &F,
        const Eigen::MatrixXi &HeOpp,
        int v, 
        int ingoing_he, 
        std::vector<std::tuple<int, int, bool>> &helist,
        int &outgoing_he)
    {
        const auto nf = F.rows();
        helist.clear();

        int f, k, i, j;
        int iHe = ingoing_he;
        do
        {
            f = iHe % nf;
            j = iHe / nf;
            k = (j + 1) % 3;
            i = (k + 1) % 3;
            if(F(f, i) != v)
            {
                return -1;
            }

            helist.emplace_back(std::make_tuple(f, j, false));
            outgoing_he = f + k * nf;
            iHe = HeOpp(f, k);
            if(iHe == ingoing_he)
            {
                return 1;
            }
        } while (iHe != -1);
        helist.emplace_back(std::make_tuple( outgoing_he % nf, outgoing_he / nf, true));
        return 0; 
    }

    DGP_INLINE bool vertex_halfedge_adjacency(
        const std::vector<std::vector<int> > &VF, 
        const std::vector<std::vector<int> > &VFi,
        const Eigen::MatrixXi &F,
        const Eigen::MatrixXi &HeOpp,
        std::vector<std::pair<bool, std::vector<std::tuple<int, int, bool>>>> &VHe)
    {
        const auto nf = F.rows();
        for(size_t v = 0; v < VF.size(); ++v)
        {
            if(VF[v].empty()) 
            {
                return false;
            }

            int f = VF[v].front(); 
            int i = VFi[v].front(); // i -> e(j, k);
            int j = (i + 1) % 3;   // j -> e(k, i);
            int k = (j + 1) % 3;   // k -> e(i, j);

            std::vector<std::tuple<int, int, bool>> helist;
            int outgoing_he; 
            const auto return_value = vertex_traversal_find_outgoing(F, HeOpp, v, f + j * nf, helist, outgoing_he);
            if(return_value == -1)
            {
                return false;
            }
            
            if(return_value == 1)
            {
                VHe.emplace_back(std::make_pair(true, helist));
                continue;
            }

            if(return_value == 0)
            {
                int _1; 
                vertex_traversal_find_ingoing(F, HeOpp, v, outgoing_he, helist, _1);
                VHe.emplace_back(std::make_pair(false, helist));
            }
        }
        return true; 
    }
}

#ifdef DGP_STATIC_LIBRARY

#endif
