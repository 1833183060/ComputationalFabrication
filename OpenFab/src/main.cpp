#include "GCodeConverter.hpp"
#include "FabSlicer.hpp"
#include "generate_commands.hpp"
#include "server.hpp"

// adds support for OBJ file reading
#include "igl/readOBJ.h"

// sample main code for FabSlicer, you should create your own main code
// int main() {
//     // load a stl file and convert to a connected mesh (obj)
//     mesh::TriMesh<double> tri_mesh(std::string(PROJECT_SOURCE_DIR) + "/data/bunny_watertight.stl", 0.01);

//     // you can visualize your mesh as an obj file to check if it's loaded correctly
//     // tri_mesh.WriteToObj(std::string(PROJECT_SOURCE_DIR) + "/data/bunny_watertight.obj");

//     // create a FabSlicer instance
//     fab_translation::FabSlicer<double> fab(tri_mesh, 0.0, 2.0, 0.03, 0.05);

//     std::vector<std::vector<std::vector<Eigen::Vector3d>>> contour;
//     std::vector<std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>>> infill_edges;

//     fab.RunTranslation(contour, infill_edges);

//     // visualize your results with ply format
//     std::string contour_file = std::string(PROJECT_SOURCE_DIR) + "/data/bunny-contour.ply";
//     fab.VisualizeContour(contour_file, 0.001, contour);

//     std::string infill_file = std::string(PROJECT_SOURCE_DIR) + "/data/bunny-infill.ply";
//     fab.VisualizeInfill(infill_file, 0.001, infill_edges);

//     return 0;
// }

// sample code for GCode generation
// int main() {
//     std::vector<std::vector<Eigen::Vector3d>> paths;
//     paths.clear();
//     // reset UI
//     network_communication::GenerateCommands::ResetCommand();
//     // translate to gcode
//     std::vector<Eigen::Vector3d> path;
//     path.clear();
//     for (int i = 0;i < 50;++i)
//         path.push_back(Eigen::Vector3d(sin(i * 0.1), cos(i * 0.1), i * 0.01));
//     paths.push_back(path);
//     fab_translation::GCodeConverter::ConvertToGCode(paths, nullptr);

//     return 0;
// }

/* Implement your code here */

// Task 0
// int main(int argc, char **argv) {
//     // load obj model
//     std::vector<std::vector<double>> verts, texCoords, norms;
//     std::vector<std::vector<int>> faces, faceTexCoords, faceNorms;
//     igl::readOBJ(std::string(PROJECT_SOURCE_DIR) + "../models/teapot.obj",
//         verts, texCoords, norms, faces, faceTexCoords, faceNorms);

//     // reset UI
//     network_communication::GenerateCommands::ResetCommand();

//     // print the cup
//     const double scale = 0.2;

//     std::vector<std::vector<Eigen::Vector3d>> paths;
//     paths.clear();
//     std::vector<Eigen::Vector3d> path;
//     path.clear();

//     for (auto &face : faces) {
//         for (int &vertIdx : face)
//             path.emplace_back(Eigen::Vector3d(verts[vertIdx][0], verts[vertIdx][1], verts[vertIdx][2]) * scale);
//         paths.push_back(path);
//         path.clear();
//     }
//     fab_translation::GCodeConverter::ConvertToGCode(paths, nullptr);

//     return 0;
// }

// Task 1-3
int main(int argc, char **argv) {
    // load input model
    mesh::TriMesh<double> tri_mesh(std::string(PROJECT_SOURCE_DIR) + "/data/bunny_watertight.stl", 0.01);
    
    // you can visualize your mesh as an obj file to check if it's loaded correctly
    // tri_mesh.WriteToObj(std::string(PROJECT_SOURCE_DIR) + "/data/bunny_watertight.obj");

    // create a FabSlicer instance
    fab_translation::FabSlicer<double> fab(tri_mesh, 0.0, 2.0, 0.03, 0.05);

    std::vector<std::vector<std::vector<Eigen::Vector3d>>> contour;
    std::vector<std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>>> infill_edges;

    fab.RunTranslation(contour, infill_edges);

    // visualize your results with ply format
    std::string contour_file = std::string(PROJECT_SOURCE_DIR) + "/data/bunny-contour.ply";
    fab.VisualizeContour(contour_file, 0.001, contour);

    std::string infill_file = std::string(PROJECT_SOURCE_DIR) + "/data/bunny-infill.ply";
    fab.VisualizeInfill(infill_file, 0.001, infill_edges);

    return 0;
}
