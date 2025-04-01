#include <dgp/smoother.h>
#include <dgp/read_labels.h>

#include <igl/opengl/glfw/imgui/ImGuiPlugin.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>

#include <igl/read_triangle_mesh.h>
#include <igl/remove_duplicate_vertices.h>

int iLabel = 0;
dgp::alg::PrecomputeParam prec_param;
dgp::alg::SDFParam sdf_param;
dgp::alg::VectorTransportParam vt_param;
double line_length = 0.5;
dgp::alg::DiffusionParam diff_param;
double phi = 0.;

Eigen::MatrixXd V; 
Eigen::MatrixXi F; 

Eigen::RowVector3d red(igl::MAYA_RED(0), igl::MAYA_RED(1), igl::MAYA_RED(2));
Eigen::RowVector3d grey(igl::MAYA_GREY(0), igl::MAYA_GREY(1), igl::MAYA_GREY(2));
Eigen::RowVector3d green(igl::MAYA_GREEN(0), igl::MAYA_GREEN(1), igl::MAYA_GREEN(2));

bool key_swith = true; 

int main(int argc, char **argv)
{
    const auto codeToLabel = [](int code, bool upper)->int
    {
        if(code == 0) return 0; 
        int pos = code % 8 == 0 ? 8 :code % 8;
        int cor = (code > 8 ? 2 : 1) * 10 + (upper ? 0 : 20);
        return pos + cor; 
    };

    const auto labelToCode = [](int label)->int
    {
        if(label == 0) return 0; 
        int pos = label % 10;
        int cor = label / 10 % 10;
        return pos + (cor == 2 || cor == 4 ? 8: 0); 
    };

    igl::opengl::glfw::Viewer viewer;
    igl::opengl::glfw::imgui::ImGuiPlugin plugin;
    viewer.plugins.push_back(&plugin);
    igl::opengl::glfw::imgui::ImGuiMenu menu;
    plugin.widgets.push_back(&menu);

    Eigen::MatrixXd inputV; 
    Eigen::MatrixXi inputF;
    std::string meshname(argv[1]);
    if(!igl::read_triangle_mesh(meshname,inputV,inputF)) 
    {
      std::cout << "Failed to load mesh " << meshname << std::endl;
      return 0;
    }
    std::string dirname, basename, extension, filename; 
    igl::pathinfo(meshname, dirname, basename, extension, filename);
    if(extension == "stl")
    {
        Eigen::MatrixXd tempV = inputV; 
        Eigen::VectorXi SVI, SVJ; 
        igl::remove_duplicate_vertices(tempV,0,inputV,SVI,SVJ);
        std::for_each(inputF.data(),inputF.data()+inputF.size(),[&SVJ](int & f){f=SVJ(f);});
    }
  
    std::string labelname(argv[2]);
    std::vector<int> vLabel;
    dgp::read_labels(labelname, vLabel);
    if(vLabel.size() != inputF.rows()) 
    {
      std::cout << "Failed to load label " << labelname << std::endl;
      return 0;
    }
    std::vector<std::string> vsLabel; 
    std::set<int> mLabel(vLabel.begin(), vLabel.end());
    std::for_each(mLabel.begin(), mLabel.end(), [&](int label) { if(label != 0) vsLabel.push_back(std::to_string(label)); });

    dgp::alg::Smoother smoother;
    smoother.load(inputV, inputF, Eigen::Map<Eigen::VectorXi>(vLabel.data(), vLabel.size()).eval());
    
    menu.callback_draw_custom_window = [&]()
    {
        {
            ImGui::SetNextWindowPos(ImVec2(180.f * menu.menu_scaling(), 10), ImGuiCond_FirstUseEver);
            ImGui::SetNextWindowSize(ImVec2(210, 50), ImGuiCond_FirstUseEver);
            ImGui::Begin(
                "Base", nullptr,
                ImGuiWindowFlags_NoSavedSettings
            );
            ImGui::Combo("label", &iLabel, vsLabel);
            ImGui::End();
        }

        {
            ImGui::SetNextWindowPos(ImVec2(180.f * menu.menu_scaling(), 10 + 50), ImGuiCond_FirstUseEver);
            ImGui::SetNextWindowSize(ImVec2(210, 130), ImGuiCond_FirstUseEver);
            ImGui::Begin(
                "precompute", nullptr,
                ImGuiWindowFlags_NoSavedSettings
            );
            ImGui::InputInt("fuzzy_kring", &prec_param.fuzzy_kring, 0, 0);
            ImGui::InputInt("block_kring", &prec_param.block_kring, 0, 0);
            ImGui::Checkbox("use_intrinsic", &prec_param.use_intrinsic);
            if(ImGui::Button("precompute"))
            {
              return smoother.precompute(std::stoi(vsLabel[iLabel]), prec_param);
            }
            ImGui::End();
        }

        {
            ImGui::SetNextWindowPos(ImVec2(180.f * menu.menu_scaling(), 10 + 50 + 130), ImGuiCond_FirstUseEver);
            ImGui::SetNextWindowSize(ImVec2(210, 180), ImGuiCond_FirstUseEver);
            ImGui::Begin(
                "sdf solve", nullptr,
                ImGuiWindowFlags_NoSavedSettings
            );
            ImGui::InputDouble("time step", &sdf_param.t, 0, 0, "%.4f");
            ImGui::Checkbox("preserve source normals", &sdf_param.preserve_source_normals);
            ImGui::Checkbox("hard constriant", &sdf_param.hard_constriant);
            ImGui::InputDouble("level set", &sdf_param.level_set, 0, 0, "%.4f");
            ImGui::InputDouble("consistency ratio", &sdf_param.consistency_ratio, 0, 0, "%.4f");
            if(ImGui::Button("sdf solve"))
            {
                return smoother.sdf_solve(std::stoi(vsLabel[iLabel]), sdf_param);
            }
            ImGui::End();
        }

        {
            ImGui::SetNextWindowPos(ImVec2(180.f * menu.menu_scaling(), 10 + 50 + 130 + 180), ImGuiCond_FirstUseEver);
            ImGui::SetNextWindowSize(ImVec2(210, 140), ImGuiCond_FirstUseEver);
            ImGui::Begin(
                "vector solve", nullptr,
                ImGuiWindowFlags_NoSavedSettings
            );
            ImGui::InputDouble("time step", &vt_param.t, 0, 0, "%.4f");
            ImGui::Checkbox("preserve source normals", &vt_param.preserve_source_normals);
            ImGui::InputInt("min spacing source normals", &vt_param.min_spacing_source_normals, 0, 0);
            ImGui::InputDouble("line length", &line_length, 0, 0, "%.4f");
            if(ImGui::Button("vector solve"))
            {
                return smoother.vector_field_solve(std::stoi(vsLabel[iLabel]), vt_param);
            }
            ImGui::End();
        }

        {
            ImGui::SetNextWindowPos(ImVec2(180.f * menu.menu_scaling(), 10 + 50 + 130 + 180 + 140), ImGuiCond_FirstUseEver);
            ImGui::SetNextWindowSize(ImVec2(210, 140), ImGuiCond_FirstUseEver);
            ImGui::Begin(
                "phi solve", nullptr,
                ImGuiWindowFlags_NoSavedSettings
            );
            ImGui::InputDouble("lambda", &diff_param.lambda, 0, 0, "%.4f");
            ImGui::InputDouble("adjacent point potential", &diff_param.adjacent_point_potential, 0, 0);
            ImGui::InputInt("min spacing constriants", &diff_param.min_spacing_constriants, 0, 0);
            ImGui::InputDouble("phi", &phi);
            if(ImGui::Button("phi solve"))
            {
                return smoother.scalar_field_solve(std::stoi(vsLabel[iLabel]), diff_param);
            }
            ImGui::End();
        }
        return true;
    };
  
    viewer.data().set_mesh(inputV,inputF);
    viewer.data().show_lines = false;

    viewer.callback_key_pressed = [&](
        igl::opengl::glfw::Viewer& viewer,
        unsigned char key,
        int modifier) -> bool
    {
        switch (key)
        {
            case '1':
            {
                Eigen::VectorXi faces_label; 
                smoother.get_mesh(std::stoi(vsLabel[iLabel]), V, F);
                smoother.get_faces_label(std::stoi(vsLabel[iLabel]), faces_label);
                Eigen::MatrixXd colors(F.rows(), 3);    
                for(int f = 0; f < F.rows(); ++f)
                {
                    colors.row(f) = static_cast<bool>(faces_label(f))? red : grey;
                }

                viewer.data().clear();
                viewer.data().set_mesh(V,F);
                viewer.core().align_camera_center(V, F);
                viewer.data().set_colors(colors);
                return true;
            }
            case '2':
            {
                Eigen::VectorXi faces_label; 
                smoother.get_mesh(std::stoi(vsLabel[iLabel]), V, F);
                smoother.get_faces_label(std::stoi(vsLabel[iLabel]), faces_label);
                viewer.data().clear();
                viewer.data().set_mesh(V,F);
                viewer.core().align_camera_center(V, F);

                Eigen::VectorXd SDF = Eigen::VectorXd::Zero(V.rows());
                smoother.get_sdf(std::stoi(vsLabel[iLabel]), SDF);
                if(key_swith)
                {
                    key_swith = !key_swith;
                    viewer.data().set_data(SDF);
                }
                else 
                {
                    key_swith = !key_swith;

                    Eigen::MatrixXd colors(F.rows(), 3);    
                    for(int f = 0; f < F.rows(); ++f)
                    {
                        colors.row(f) = static_cast<bool>(faces_label(f))? red : grey;
                    }
                    viewer.data().set_colors(colors);

                    Eigen::VectorXd level_sets;
                    level_sets.resize(1);
                    level_sets << sdf_param.level_set;
                    Eigen::VectorXi iI;
                    Eigen::MatrixXd iV;
                    Eigen::MatrixXi iE; 
                    igl::isolines(V, F, SDF, level_sets, iV, iE, iI);
                    viewer.data().set_edges(iV, iE, green);
                }
                return true;
            }
            case '3':
            {

                Eigen::VectorXi faces_label; 
                smoother.get_mesh(std::stoi(vsLabel[iLabel]), V, F);
                smoother.get_faces_label(std::stoi(vsLabel[iLabel]), faces_label);
                Eigen::MatrixXd X; 
                smoother.get_flow(std::stoi(vsLabel[iLabel]), X);
                Eigen::MatrixXd colors(F.rows(), 3);    
                for(int f = 0; f < F.rows(); ++f)
                {
                    colors.row(f) = static_cast<bool>(faces_label(f))? red : grey;
                }

                viewer.data().clear();
                viewer.data().set_mesh(V,F);
                viewer.core().align_camera_center(V, F);
                viewer.data().set_colors(colors);
                viewer.data().add_edges(V, V + line_length * X.rowwise().normalized(), green);
                return true;
            }
            case '4':
            {
                Eigen::VectorXi faces_label; 
                smoother.get_mesh(std::stoi(vsLabel[iLabel]), V, F);
                smoother.get_faces_label(std::stoi(vsLabel[iLabel]), faces_label);
                viewer.data().clear();
                viewer.data().set_mesh(V,F);
                viewer.core().align_camera_center(V, F);

                Eigen::VectorXd PHI = Eigen::VectorXd::Zero(V.rows());
                smoother.get_phi(std::stoi(vsLabel[iLabel]), PHI);
                if(key_swith)
                {
                    key_swith = !key_swith;
                    viewer.data().set_data(PHI);
                }
                else 
                {
                    key_swith = !key_swith;

                    Eigen::MatrixXd colors(F.rows(), 3);    
                    for(int f = 0; f < F.rows(); ++f)
                    {
                        colors.row(f) = static_cast<bool>(faces_label(f))? red : grey;
                    }
                    viewer.data().set_colors(colors);

                    Eigen::VectorXd level_sets;
                    level_sets.resize(1);
                    level_sets << phi;
                    Eigen::VectorXi iI;
                    Eigen::MatrixXd iV;
                    Eigen::MatrixXi iE; 
                    igl::isolines(V, F, PHI, level_sets, iV, iE, iI);
                    viewer.data().set_edges(iV, iE, green);
                }
                return true;
            }
            default: return false; 
        }
        return true; 
    };
    viewer.launch();
    return 0;
}