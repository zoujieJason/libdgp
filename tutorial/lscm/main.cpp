#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <igl/boundary_loop.h>
#include <igl/opengl/glfw/Viewer.h>

#include <dgp/lscm.h>

int main(int argc, char **argv)
{
    Eigen::MatrixXd V; 
    Eigen::MatrixXi F; 
    if(!igl::readOBJ(LIBDGP_DATA_PATH"/head.obj", V, F))
    {
        return 0;
    }

    Eigen::MatrixXd UV; 
    dgp::lscm(V, F, UV);
    UV.conservativeResize(UV.rows(), 3);
    UV.col(2).setZero();

    igl::writeOBJ(LIBDGP_DATA_PATH"/head_lscm.obj", UV, F);

    igl::opengl::glfw::Viewer viewer;
    viewer.data().set_mesh(UV, F);
    viewer.data().set_face_based(true);
    viewer.launch();
    return 0;
}