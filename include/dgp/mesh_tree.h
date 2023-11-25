#ifndef LIBDGP_MESH_TREE_H
#define LIBDGP_MESH_TREE_H

#include "dgp_inline.h"

#include <Eigen/Dense>

#include <vector>

namespace dgp
{
    class MeshTree
    {
    public:
        struct Node
        {
            Node(int fid = -1, int cid = -1) : fid_(fid), cid_(cid) {}
            int fid_;
            int cid_;
            std::weak_ptr<Node> parent_;
            std::vector<std::weak_ptr<Node> > children_;
        };

    public:
        MeshTree() = default;
        void build(int nv, const Eigen::MatrixXi &F);
        std::vector<std::weak_ptr<Node>>::iterator begin();
        std::vector<std::weak_ptr<Node>>::iterator end();
        void traverse(std::weak_ptr<Node> node, const std::function<void(std::weak_ptr<Node>)> &apply) const;

    private:
        std::vector<std::weak_ptr<Node>> heads_;
        std::vector<std::shared_ptr<Node>> nodes_;   
    };
}

#ifndef DGP_STATIC_LIBRARY
#include   "mesh_tree.cpp"
#endif

#endif //LIBDGP_MESH_TREE_H