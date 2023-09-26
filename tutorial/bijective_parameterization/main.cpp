#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <igl/boundary_loop.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/doublearea.h>

#include <dgp/unit_circle_points.h>
#include <dgp/barycentric_embedding.h>
#include <dgp/bijective_parameterization.h>

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
    dgp::unit_circle_points(mesh_boundary.size(), 10.0, embedding_boundary);

    Eigen::MatrixXd UV; 
    if(!dgp::barycentric_embedding(V, F, mesh_boundary, embedding_boundary, UV))
    {
        return 0;
    }
    
    dgp::BijectiveParameterization BP; 
    BP.initialize(V, F, UV);
    
    int k = 0; 
    double g_norm, lambda; 
    BP.per_iteration(k, g_norm, lambda);
    
    igl::opengl::glfw::Viewer viewer;
    viewer.data().set_mesh(UV, F);
    viewer.data().set_face_based(true);
    viewer.launch();
    return 0;
}