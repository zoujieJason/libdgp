#include <dgp/smoother.h>
#include <dgp/read_labels.h>

#include <igl/opengl/glfw/imgui/ImGuiPlugin.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/opengl/glfw/imgui/SelectionWidget.h>

#include <igl/read_triangle_mesh.h>
#include <igl/remove_duplicate_vertices.h>
#include <igl/AABB.h>
#include <igl/screen_space_selection.h>

int iLabel = 0;
dgp::alg::PrecomputeParam prec_param;
dgp::alg::SDFParam sdf_param;
dgp::alg::VectorTransportParam vt_param;
double line_length = 0.5;
dgp::alg::DiffusionParam diff_param;
double phi = 0.;

Eigen::MatrixXd V; 
Eigen::MatrixXi F; 
igl::AABB<Eigen::MatrixXd, 3> tree;
Eigen::VectorXd W; 
Eigen::Array<double,Eigen::Dynamic,1> and_visible;
Eigen::VectorXd S; 
std::vector<size_t> selected_indices; 
double weight = 1.;
double potential = 0.01;

Eigen::RowVector3d red(igl::MAYA_RED(0), igl::MAYA_RED(1), igl::MAYA_RED(2));
Eigen::RowVector3d grey(igl::MAYA_GREY(0), igl::MAYA_GREY(1), igl::MAYA_GREY(2));
Eigen::RowVector3d green(igl::MAYA_GREEN(0), igl::MAYA_GREEN(1), igl::MAYA_GREEN(2));
bool key_swith = true; 

int main(int argc, char **argv)
{
    igl::opengl::glfw::Viewer viewer;
    igl::opengl::glfw::imgui::ImGuiPlugin plugin;
    viewer.plugins.push_back(&plugin);
    igl::opengl::glfw::imgui::ImGuiMenu menu;
    plugin.widgets.push_back(&menu);
    igl::opengl::glfw::imgui::SelectionWidget widget;
    plugin.widgets.push_back(&widget);
    widget.mode = igl::opengl::glfw::imgui::SelectionWidget::Mode::OFF; 

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

    std::vector<std::vector<int> > inputVF, inputVFi; 
    igl::vertex_triangle_adjacency(inputV, inputF, inputVF, inputVFi);

    Eigen::VectorXi input_faces_label = Eigen::Map<Eigen::VectorXi>(vLabel.data(), vLabel.size());
    dgp::alg::Smoother smoother;
    smoother.load(inputV, inputF, input_faces_label);

    menu.callback_draw_custom_window = [&]()
    {
        // base
        float window_current_height = 10.f; 
        {
            float window_size_y = 80.f;
            ImGui::SetNextWindowPos(ImVec2(180.f * menu.menu_scaling(), window_current_height), ImGuiCond_FirstUseEver);
            ImGui::SetNextWindowSize(ImVec2(210, window_size_y), ImGuiCond_FirstUseEver);
            ImGui::Begin("Base", nullptr, ImGuiWindowFlags_NoSavedSettings);
            ImGui::Combo("label", &iLabel, vsLabel);
            if(ImGui::Button("clear"))
            {
                selected_indices.clear();
                prec_param.remove_vertices.clear();
                prec_param.block_vertices.clear();
                vt_param.block_vertices.clear();
                diff_param.block_vertices.clear();
                diff_param.potential_vertices.clear();
                smoother.load(inputV, inputF, Eigen::Map<Eigen::VectorXi>(vLabel.data(), vLabel.size()).eval());
            }
            ImGui::End();
            window_current_height += window_size_y; 
        }

        // precompute
        {
            float window_size_y = 170.f;
            ImGui::SetNextWindowPos(ImVec2(180.f * menu.menu_scaling(), window_current_height), ImGuiCond_FirstUseEver);
            ImGui::SetNextWindowSize(ImVec2(210, window_size_y), ImGuiCond_FirstUseEver);
            ImGui::Begin("precompute", nullptr, ImGuiWindowFlags_NoSavedSettings);
            ImGui::InputInt("fuzzy_kring", &prec_param.fuzzy_kring, 0, 0);
            ImGui::InputInt("block_kring", &prec_param.block_kring, 0, 0);
            ImGui::Checkbox("use_intrinsic", &prec_param.use_intrinsic);
            if(ImGui::Button("precompute"))
            {
                smoother.precompute(std::stoi(vsLabel[iLabel]), prec_param);
                smoother.get_mesh(std::stoi(vsLabel[iLabel]), V, F);
                and_visible.setZero(V.rows());
                W.setZero(V.rows());
                S.setZero(V.rows());
                tree.init(V, F);
                prec_param.remove_vertices.clear();
            }
            if(ImGui::Button("save remove vertices"))
            {
                prec_param.remove_vertices = selected_indices; 
                std::cout << "selected " << prec_param.remove_vertices.size() << " remove vertices\n";
            }
            if(ImGui::Button("save block vertices"))
            {
                prec_param.block_vertices = selected_indices; 
                std::cout << "selected " << prec_param.block_vertices.size() << " block vertices\n";
            }
            ImGui::End();
            window_current_height += window_size_y; 
        }

        // sdf solve
        {
            float window_size_y = 170.f;
            ImGui::SetNextWindowPos(ImVec2(180.f * menu.menu_scaling(), window_current_height), ImGuiCond_FirstUseEver);
            ImGui::SetNextWindowSize(ImVec2(210, window_size_y), ImGuiCond_FirstUseEver);
            ImGui::Begin("sdf solve", nullptr, ImGuiWindowFlags_NoSavedSettings);
            ImGui::InputDouble("time step", &sdf_param.t, 0, 0, "%.4f");
            ImGui::Checkbox("preserve source normals", &sdf_param.preserve_source_normals);
            ImGui::Checkbox("hard constriant", &sdf_param.hard_constriant);
            ImGui::InputDouble("level set", &sdf_param.level_set, 0, 0, "%.4f");
            ImGui::InputDouble("consistency ratio", &sdf_param.consistency_ratio, 0, 0, "%.4f");
            if(ImGui::Button("sdf solve"))
            {
                smoother.sdf_solve(std::stoi(vsLabel[iLabel]), sdf_param); 
            }
            ImGui::End();
            window_current_height += window_size_y; 
        }

        // vector solve
        {
            float window_size_y = 170.f;
            ImGui::SetNextWindowPos(ImVec2(180.f * menu.menu_scaling(), window_current_height), ImGuiCond_FirstUseEver);
            ImGui::SetNextWindowSize(ImVec2(210, window_size_y), ImGuiCond_FirstUseEver);
            ImGui::Begin("vector solve", nullptr, ImGuiWindowFlags_NoSavedSettings);
            ImGui::InputDouble("time step", &vt_param.t, 0, 0, "%.4f");
            ImGui::Checkbox("preserve source normals", &vt_param.preserve_source_normals);
            ImGui::InputInt("min spacing source normals", &vt_param.min_spacing_source_normals, 0, 0);
            ImGui::InputDouble("line length", &line_length, 0, 0, "%.4f");
            if(ImGui::Button("vector solve"))
            {
                smoother.vector_field_solve(std::stoi(vsLabel[iLabel]), vt_param); 
            }
            if(ImGui::Button("save block vertices"))
            {
                vt_param.block_vertices = selected_indices; 
                std::cout << "vt_param: select " << vt_param.block_vertices.size() << " as block vertices\n";
            }
            ImGui::End();
            window_current_height += window_size_y; 
        }

        // phi solve
        {
            float window_size_y = 260.f;
            ImGui::SetNextWindowPos(ImVec2(180.f * menu.menu_scaling(), window_current_height), ImGuiCond_FirstUseEver);
            ImGui::SetNextWindowSize(ImVec2(210, window_size_y), ImGuiCond_FirstUseEver);
            ImGui::Begin("phi solve", nullptr, ImGuiWindowFlags_NoSavedSettings);
            ImGui::InputDouble("lambda", &diff_param.lambda, 0, 0, "%.4f");
            ImGui::InputDouble("adjacent point potential", &diff_param.adjacent_point_potential, 0, 0);
            ImGui::InputInt("min spacing constriants", &diff_param.min_spacing_constriants, 0, 0);
            ImGui::InputDouble("phi", &phi);
            if(ImGui::Button("phi solve"))
            {
                smoother.scalar_field_solve(std::stoi(vsLabel[iLabel]), diff_param); 
            }
            if(ImGui::Button("save block vertices"))
            {
                diff_param.block_vertices.clear();
                diff_param.block_vertices.reserve(selected_indices.size());
                for(auto selected_index: selected_indices)
                {
                    diff_param.block_vertices.emplace_back(std::make_pair(selected_index, weight));
                }
                std::cout << "diff_param: select " << diff_param.block_vertices.size() << " as block vertices\n";
            }
            ImGui::InputDouble("weight", &weight);
            ImGui::InputDouble("lambda regularization", &diff_param.lambda_regularization, 0, 0, "%.4f");
            ImGui::InputDouble("potential", &potential);
            if(ImGui::Button("save potential vertices"))
            {
                diff_param.potential_vertices.clear();
                diff_param.potential_vertices.reserve(selected_indices.size());
                for(auto selected_index: selected_indices)
                {
                    diff_param.potential_vertices.emplace_back(std::make_pair(selected_index, potential));
                }
                std::cout << "select " << diff_param.potential_vertices.size() << " as potential vertices\n";
            }
            ImGui::End();
            window_current_height += window_size_y; 
        }
    };

    widget.callback = [&]()
    {
        igl::screen_space_selection(
            V,F,tree,
            viewer.core().view,
            viewer.core().proj,
            viewer.core().viewport,
            widget.L,
            W,
            and_visible);

        S = W;
        S.array() *= and_visible;
        selected_indices.clear();
        for(Eigen::Index i = 0; i < S.size(); ++i)
        {
            if(S(i)>0.)
            {
                selected_indices.push_back(i); 
            }
        }
        W.setZero(V.rows());
        and_visible.setZero(V.rows());
    };

    viewer.data().set_mesh(inputV,inputF);
    viewer.data().show_lines = false;
    viewer.core().set_rotation_type(igl::opengl::ViewerCore::ROTATION_TYPE_TRACKBALL);

    viewer.callback_key_pressed = [&](
        igl::opengl::glfw::Viewer& viewer,
        unsigned char key,
        int modifier) -> bool
    {
        switch (key)
        {
            case '0':
            {
                const Eigen::MatrixXd CM = (Eigen::MatrixXd(2,3)<<
                0.3,0.3,0.5,                 
                255.0/255.0,228.0/255.0,58.0/255.0).finished();

                viewer.data().clear();
                viewer.data().set_mesh(V,F);
                viewer.core().align_camera_center(V, F);
                viewer.data().set_data(S,0,1,igl::COLOR_MAP_TYPE_PLASMA,2);
                viewer.data().set_colormap(CM);
                viewer.data().set_points(V(selected_indices, Eigen::all), red);
                return true;
            }
            case '1':
            {
                Eigen::VectorXi faces_label = Eigen::VectorXi::Zero(F.rows());
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