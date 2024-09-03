#ifndef LIBDGP_EXTRACT_KRING_H
#define LIBDGP_EXTRACT_KRING_H

#include "dgp_inline.h"

#include <Eigen/Dense>

#include <vector>
#include <map>

namespace dgp
{
    template<
        typename DerivedF,
        typename VFType,
        typename VType,
        typename KType, 
        typename KRingType> 
    DGP_INLINE void extract_kring(
        const Eigen::PlainObjectBase<DerivedF> & F,
        const std::vector<std::vector<VFType>> & VF, 
        VType v,
        KType k,
        std::vector<KRingType> & kring_vertices,
        std::vector<KRingType> & kring_faces)
    {
        std::map<VType, int> D; 
        kring_vertices.push_back(v); 
        D[v] = 0;

        VType current_probe = 0; 
        KType level;
        while(current_probe < kring_vertices.size() && (level = D[kring_vertices[current_probe]]) < k)
        {
            v = kring_vertices[current_probe++];
            for(auto vf: VF[v])
            {
                for(int i = 0; i < 3; ++i)
                {
                    if(D.insert(std::make_pair(F(vf, i), level + 1)).second)
                    {
                        kring_vertices.push_back(F(vf, i));
                    }
                }
                kring_faces.push_back(vf);
            }
        }
    }
}

#endif //LIBDGP_EXTRACT_KRING_H