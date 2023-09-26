#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <igl/boundary_loop.h>
#include <igl/opengl/glfw/Viewer.h>

#include <dgp/unit_circle_points.h>
#include <dgp/barycentric_embedding.h>

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
    UV.conservativeResize(UV.rows(), 3);
    UV.col(2).setZero();

    igl::writeOBJ(LIBDGP_DATA_PATH"/head_uv.obj", UV, F);

    igl::opengl::glfw::Viewer viewer;
    viewer.data().set_mesh(UV, F);
    viewer.data().set_face_based(true);
    viewer.launch();
    return 0;
}