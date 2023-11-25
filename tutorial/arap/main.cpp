#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <igl/boundary_loop.h>
#include <igl/opengl/glfw/Viewer.h>

#include <dgp/unit_circle_points.h>
#include <dgp/barycentric_embedding.h>

#include <dgp/arap.h>

int main(int argc, char **argv)
{
    Eigen::MatrixXd V; 
    Eigen::MatrixXi F; 
    if(!igl::readOBJ(LIBDGP_DATA_PATH"/head.obj", V, F))
    {
        return 0;
    }

    std::vector<size_t> mesh_boundary; 
    igl::boundary_loop(F, mesh_boundary);

    Eigen::MatrixXd embedding_boundary; 
    dgp::unit_circle_points(mesh_boundary.size(), 1.0, embedding_boundary);

    Eigen::MatrixXd UV; 
    if(!dgp::barycentric_embedding(V, F, mesh_boundary, embedding_boundary, UV))
    {
        return 0;
    }

    Eigen::MatrixXd UV_arap;
    dgp::arap(V, F, UV, UV_arap);

    igl::opengl::glfw::Viewer viewer;
    viewer.data().set_mesh(UV_arap, F);
    viewer.data().set_face_based(true);
    viewer.launch();
    return 0;
}