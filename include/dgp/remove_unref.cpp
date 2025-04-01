#include "remove_unref.h"

#include <igl/remove_unreferenced.h>

namespace dgp
{
    void remove_unref(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& F,
        const Eigen::VectorXi& K,
        Eigen::MatrixXd& Vo,
        Eigen::MatrixXi& Fo,
        Eigen::VectorXi& Io,
        Eigen::VectorXi& Jo,
        Eigen::VectorXi& FoMAP)
    {
        auto extract_faces = (K.array() == 1).count();
        FoMAP.setZero(extract_faces);
        Eigen::MatrixXi Fi(extract_faces, 3);
        for (int i = 0, f = 0; i < F.rows(); i++)
        {
            if (K(i) == 1)
            {
                FoMAP(f) = i;
                Fi.row(f++) = F.row(i);
            }
        }

        igl::remove_unreferenced(V, Fi, Vo, Fo, Io, Jo);
    }
}