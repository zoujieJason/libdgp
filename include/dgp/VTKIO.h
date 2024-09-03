#ifndef LIBDGP_VTKIO_H
#define LIBDGP_VTKIO_H

#include "dgp_inline.h"

#include <Eigen/Dense>

#include <string>

namespace dgp
{
    enum class VTK_DATASET
    {
        STRUCTURED_POINTS,
        STRUCTURED_GRID,
        UNSTRUCTURED_GRID,
        POLYDATA,
        RECTILINEAR_GRID,
        FIELD
    };

    enum class VTK_POLYDATA_TYPE
    {
        POINTS,
        VERTICES,
        LINES,
        POLYGONS,
        TRIANGLE_STRIPS
    };

    enum class VTK_ATTRIBUTE_TYPE
    {
        POINT_DATA,
        CELL_DATA
    };

    enum class VTK_ATTRIBUTEDATA_TYPE
    {
        SCALARS,
        VECTORS,
        NORMALS,
        TEXTURE_COORDINATES,
        TENSORS,
        //FIELD
    };

    class VTKIO
    {
    public:
        template<typename OS>
        DGP_INLINE static void InitHead(
            OS &os, 
            VTK_DATASET data_set, 
            std::string title = "legacy format",
            std::string version = "2.0",
            std::string type = "ASCII")
        {
            os << "# vtk DataFile Version " << version << std::endl;
            os << title << std::endl;
            os << type << std::endl;
            switch (data_set)
            {
                case VTK_DATASET::STRUCTURED_POINTS:
                    os << "DATASET STRUCTURED_POINTS" << std::endl;
                    break;
                case VTK_DATASET::STRUCTURED_GRID:
                    os << "DATASET STRUCTURED_GRID" << std::endl;
                    break;
                case VTK_DATASET::UNSTRUCTURED_GRID:
                    os << "DATASET UNSTRUCTURED_GRID" << std::endl;
                    break;
                case VTK_DATASET::POLYDATA:
                    os << "DATASET POLYDATA" << std::endl;
                    break;
                case VTK_DATASET::RECTILINEAR_GRID:  
                    os << "DATASET RECTILINEAR_GRID" << std::endl;
                    break;          
                case VTK_DATASET::FIELD:     
                    os << "DATASET FIELD" << std::endl;
                    break;
            }
        }

        template<typename OS, typename Index>
        DGP_INLINE static void InitAttribute(
            VTK_ATTRIBUTE_TYPE type,
            OS &os,
            Index n)
        {
            switch (type)
            {
                case VTK_ATTRIBUTE_TYPE::POINT_DATA:
                    os << "POINT_DATA " << n << std::endl;
                    break;
                case VTK_ATTRIBUTE_TYPE::CELL_DATA:
                    os << "CELL_DATA " << n << std::endl;
                    break;
            }
        }

        template<typename DerivedPD, typename OS>
        DGP_INLINE static void PolyData(
            const Eigen::PlainObjectBase<DerivedPD> &PD,
            VTK_POLYDATA_TYPE type,
            OS &os)
        {
            switch (type)
            {
                case VTK_POLYDATA_TYPE::POINTS: os << "POINTS "; break;
                case VTK_POLYDATA_TYPE::VERTICES: os << "VERTICES "; break;
                case VTK_POLYDATA_TYPE::LINES: os << "LINES "; break;
                case VTK_POLYDATA_TYPE::POLYGONS: os << "POLYGONS "; break;
                case VTK_POLYDATA_TYPE::TRIANGLE_STRIPS: os << "TRIANGLE_STRIPS "; break;
            }
            os << PD.rows() << " " << (type == VTK_POLYDATA_TYPE::POINTS ? "float" : std::to_string(PD.rows() + PD.size())) << std::endl;
            std::string _row_prefix = type == VTK_POLYDATA_TYPE::POINTS ? "" : std::to_string(PD.cols()) + " ";
            os << PD.format(Eigen::IOFormat(Eigen::FullPrecision, Eigen::DontAlignCols, " ", "\n", _row_prefix, "", "", "\n"));
        }

        template<typename DerivedAD, typename OS>
        DGP_INLINE static void AttributeData(
            const Eigen::PlainObjectBase<DerivedAD> &AD,
            VTK_ATTRIBUTEDATA_TYPE type,
            OS &os,
            std::string attribute_name,
            std::string lookup_table = "default")
        {
            switch (type)
            {
                case VTK_ATTRIBUTEDATA_TYPE::SCALARS: os << "SCALARS " << attribute_name << " float\nLOOKUP_TABLE " << lookup_table << "\n"; break;
                case VTK_ATTRIBUTEDATA_TYPE::VECTORS: os << "VECTORS " << attribute_name << " float\n"; break;
                case VTK_ATTRIBUTEDATA_TYPE::NORMALS: os << "NORMALS " << attribute_name << " float\n"; break;
                case VTK_ATTRIBUTEDATA_TYPE::TEXTURE_COORDINATES: os << "TEXTURE_COORDINATES " << attribute_name << " " << AD.cols() << " float\n"; break;
                case VTK_ATTRIBUTEDATA_TYPE::TENSORS: os << "TENSORS " << attribute_name << " float\n"; break;
            }
            os << AD.format(Eigen::IOFormat(Eigen::FullPrecision, Eigen::DontAlignCols, " ", "\n", "", "", "", "\n"));
        }
    };
}

#endif //LIBDGP_VTKIO_H