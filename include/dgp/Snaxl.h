#ifndef LIBDGP_SNAXL_H
#define LIBDGP_SNAXL_H

#include "dgp_inline.h"

#include <Eigen/Dense>

#include <new>

namespace dgp
{
    class Snaxl 
    {
    public:
        enum class Type { Vertex, Edge };

        struct VertexSnaxl
        {
            int v{-1};
            VertexSnaxl(int v_) : v(v_) {}
        };

        struct EdgeSnaxl
        {
            int ue{-1};
            double t{-1.}; 
            int v0{-1};
            int v1{-1};
            EdgeSnaxl(int ue_, double t_, int v0_, int v1_): ue(ue_), t(t_), v0(v0_), v1(v1_) {};
        };

        Type GetSnaxlType() const { return type_; }

        DGP_INLINE Snaxl(const VertexSnaxl &snaxl): type_(Type::Vertex) 
        {
            new (&vertex_snaxl_) VertexSnaxl(snaxl);
        }

        DGP_INLINE Snaxl(const EdgeSnaxl &snaxl): type_(Type::Edge)
        {
            new (&edge_snaxl_) EdgeSnaxl(snaxl);
        }

        Snaxl(const Snaxl& other) : type_(other.type_) 
        {
            if (type_ == Type::Vertex) 
            {
                new (&vertex_snaxl_) VertexSnaxl(other.vertex_snaxl_);
            } 
            else 
            {
                new (&edge_snaxl_) EdgeSnaxl(other.edge_snaxl_);
            }
        }

        Snaxl& operator=(const Snaxl& other) 
        {
            if (this != &other)
            {
                if (type_ == Type::Vertex) 
                {
                    vertex_snaxl_.~VertexSnaxl();
                } 
                else 
                {
                    edge_snaxl_.~EdgeSnaxl();
                }
                
                type_ = other.type_;
                
                if (type_ == Type::Vertex) 
                {
                    new (&vertex_snaxl_) VertexSnaxl(other.vertex_snaxl_);
                } 
                else 
                {
                    new (&edge_snaxl_) EdgeSnaxl(other.edge_snaxl_);
                }
            }
            return *this;
        }

        ~Snaxl() {
            if (type_ == Type::Vertex) 
            {
                vertex_snaxl_.~VertexSnaxl();
            } 
            else 
            {
                edge_snaxl_.~EdgeSnaxl();
            }
        }

        const VertexSnaxl & AsVertexSnaxl() const 
        {
            if(type_ != Type::Vertex) 
            {
                throw std::bad_cast();
            }
            return vertex_snaxl_; 
        }

        const EdgeSnaxl & AsEdgeSnaxl() const 
        {
            if(type_ != Type::Edge) 
            {
                throw std::bad_cast();
            }
            return edge_snaxl_; 
        }
        
    private:
        Type type_; 
        union 
        {
            VertexSnaxl vertex_snaxl_; 
            EdgeSnaxl edge_snaxl_;
        };
    };
}

#endif //LIBDGP_SNAXL_H