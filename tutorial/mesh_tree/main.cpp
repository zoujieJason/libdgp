#include <igl/readOBJ.h>

#include <dgp/mesh_tree.h>

int main(int argc, char **argv)
{
    Eigen::MatrixXd V; 
    Eigen::MatrixXi F; 
    if(!igl::readOBJ(LIBDGP_DATA_PATH"/head.obj", V, F))
    {
        return 0;
    }

    dgp::MeshTree mesh_tree;
    mesh_tree.build(V.rows(), F);
    for(auto head : mesh_tree)
    {
        mesh_tree.traverse(head, [](std::weak_ptr<dgp::MeshTree::Node> node)
        {
            std::cout << node.lock()->fid_ << " ";
        });
    }
    return 0;
}