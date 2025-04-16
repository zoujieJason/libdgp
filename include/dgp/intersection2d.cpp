#include "intersection2d.h"

namespace dgp
{
    DGP_INLINE bool intersection2d(
        const Eigen::Vector3d& O1, 
        const Eigen::Vector3d& E1,
        const Eigen::Vector3d& O2,
        const Eigen::Vector3d& E2,
        Eigen::Vector3d& pt, 
        double& t1, 
        double& t2)
    {
        Eigen::Vector3d D1 = E1 - O1;
        Eigen::Vector3d D2 = E2 - O2;
        Eigen::Vector3d O12 = O2 - O1;
    
        double D = -D1.x() * D2.y() + D1.y() * D2.x();
        if (std::abs(D) < 1e-8)
        {
             return false; 
        }
    
        t1 = (-O12.x() * D2.y() + O12.y() * D2.x()) / D;
        t2 = (D1.x() * O12.y() - D1.y() * O12.x()) / D;
    
        if (t1 < 0 || t1 > 1 || t2 < 0 || t2 > 1) 
        {
            return false;
        }
    
        pt = O1 + t1 * D1;
        pt.z() = 0.; 
        return true;
    }
}