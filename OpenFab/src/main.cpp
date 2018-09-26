#include "GCodeConverter.hpp"
#include "FabSlicer.hpp"
#include "generate_commands.hpp"
#include "server.hpp"

// adds support for OBJ file reading
#include "igl/readOBJ.h"

// Assignment 2
#include "voxelizer.hpp"
#include "marching_cube.hpp"
#include <iostream>
#include <ctime>
#include <string>

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

// Assignment 1, Task 0
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

// Assignment 1, Task 1-3
// int main(int argc, char **argv) {
//     // load input model
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

// Assignment 2
int main(int argc, char* argv[]) {
    // Usage:
    // Slow version:
    // ./OpenFab bunny.stl 10.0
    // ./OpenFab fandisk.stl 0.5
    // ./OpenFab spot.stl 0.25
    // ./OpenFab dragon.stl 1.0
    // Fast version:
    // ./OpenFab bunny.stl 2.0 fast
    // ./OpenFab fandisk.stl 0.05 fast
    // ./OpenFab spot.stl 0.125 fast
    // ./OpenFab dragon.stl 0.05 fast
    // Approximation:
    // ./OpenFab bunny_with_hole.stl 2.0 approx
    // ./OpenFab spot_with_whole.stl 0.125 approx
    // Marching cube version:
    // ./OpenFab bunny_voxel_info.txt
    // ./OpenFab dragon_voxel_info.txt
    // ./OpenFab fandisk_voxel_info.txt
    // ./OpenFab spot_voxel_info.txt
    if (argc == 2) {
        // Marching cube version.
        const std::string info_file(argv[1]);
        mesh::MarchingCube<double> mc(std::string(PROJECT_SOURCE_DIR) + "/data/assignment2/" + info_file);
        mc.BuildMesh();
        const std::string name = info_file.substr(0, info_file.size() - std::string("_voxel_info.txt").size());
        mc.ExportMeshToFile(name + "_mc.stl");
        return 0;
    }

    int t0 = std::clock();
    const std::string stl_name(argv[1]);
    const double dx = std::stod(argv[2]);
    mesh::Voxelizer<double> voxelizer(std::string(PROJECT_SOURCE_DIR) + "/data/assignment2/" + stl_name, dx);
    int t1 = std::clock();
    std::cout << "load mesh success... " << (double)(t1 - t0) / 1000000.0 << " seconds." << std::endl;
    std::cout << "Bounding box: " << voxelizer.pmin().transpose() << ", " << voxelizer.pmax().transpose() << std::endl;
    std::cout << "Number of voxels: " << voxelizer.voxel_num().transpose() << std::endl;
    if (argc == 3) {
        voxelizer.BasicVoxelization();
    } else {
        const std::string flag(argv[3]);
        if (flag == "fast")
            voxelizer.AdvancedVoxelization();
        else if (flag == "approx")
            voxelizer.AdvancedVoxelizationWithApproximation();
        else {
            std::cout << "ERROR: unexpected flag" << std::endl;
            exit(0);
        }
    }
    std::cout << "Voxelization done..." << std::endl;
    // Export results to mesh.
    const std::string stl_prefix = stl_name.substr(0, stl_name.size() - 4);
    const std::string voxel_file_name = stl_prefix + "_voxel.stl";
    std::cout << "Saving results to " << voxel_file_name << std::endl;
    voxelizer.WriteVoxelToMesh(voxel_file_name);
    std::cout << "Results saved..." << std::endl;
    return 0;
}
