#include "linear_material.hpp"
#include "deformable_body.hpp"
#include <stdio.h>
#include <iostream>
#include "poly_mesh.hpp"
#include "hexahedral_mesh.hpp"
#include "hex_deformable_body.hpp"
#include "typedefs.hpp"
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCholesky>

#include "voxelizer.hpp"
#include "marching_cube.hpp"

#include <stdio.h>
#include <iostream>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCholesky>
#include <fstream>
#include <set>

#include "GeometryExploration.hpp"

void write_facet(std::ofstream& file, Eigen::Vector3d p0, Eigen::Vector3d p1, Eigen::Vector3d p2) {
    file << "facet normal 0 0 0\n";
    file << "outer loop\n";
    file << "vertex " + std::to_string(p0.x()) + " " + std::to_string(p0.y()) + " " + std::to_string(p0.z()) + "\n";
    file << "vertex " + std::to_string(p1.x()) + " " + std::to_string(p1.y()) + " " + std::to_string(p1.z()) + "\n";
    file << "vertex " + std::to_string(p2.x()) + " " + std::to_string(p2.y()) + " " + std::to_string(p2.z()) + "\n";
    file << "endloop\n";
    file << "endfacet\n";
}

void write_voxel_grid(const std::string& filename, const Eigen::Matrix<double, 3, Eigen::Dynamic>& vertices, const Eigen::Matrix<int, 8, Eigen::Dynamic>& elements) {

    std::ofstream file;
    file.open(filename, std::ios::out);

    file << "solid voxel_grid \n";

    for (size_t i = 0; i < elements.cols(); ++i) {

        write_facet(file, vertices.col(elements(0, i)), vertices.col(elements(1, i)), vertices.col(elements(2, i)));
        write_facet(file, vertices.col(elements(1, i)), vertices.col(elements(3, i)), vertices.col(elements(2, i)));



        write_facet(file, vertices.col(elements(4, i)), vertices.col(elements(5, i)), vertices.col(elements(0, i)));
        write_facet(file, vertices.col(elements(5, i)), vertices.col(elements(1, i)), vertices.col(elements(0, i)));


        write_facet(file, vertices.col(elements(2, i)), vertices.col(elements(3, i)), vertices.col(elements(6, i)));
        write_facet(file, vertices.col(elements(3, i)), vertices.col(elements(7, i)), vertices.col(elements(6, i)));


        write_facet(file, vertices.col(elements(6, i)), vertices.col(elements(7, i)), vertices.col(elements(4, i)));
        write_facet(file, vertices.col(elements(7, i)), vertices.col(elements(5, i)), vertices.col(elements(4, i)));

        write_facet(file, vertices.col(elements(1, i)), vertices.col(elements(3, i)), vertices.col(elements(5, i)));
        write_facet(file, vertices.col(elements(3, i)), vertices.col(elements(7, i)), vertices.col(elements(5, i)));



        write_facet(file, vertices.col(elements(0, i)), vertices.col(elements(2, i)), vertices.col(elements(4, i)));
        write_facet(file, vertices.col(elements(2, i)), vertices.col(elements(6, i)), vertices.col(elements(4, i)));

    }
    file << "endsolid";

    file.close();

}

void SolvePerformance(std::string stl_name, double& compliance, double& num_voxels) {
    // Q4: implement the pipeline from design space to the performance space
    // push_back your results to the compliance and num_voxels
    const int dim = 3;
    materials::LinearElasticityMaterial<dim, double> linear_elasticity_material(10000000, 0.45);
    double dx = 0.25;

    // Your code starts here!!
    mesh::Voxelizer<double> voxelizer(stl_name, dx);
    voxelizer.AdvancedVoxelization();
    materials::HexahedralMesh<double> hex_mesh = std::move(voxelizer.ConvertToHexMesh());
    const auto &vertex = hex_mesh.vertex();
    const auto &element = hex_mesh.element();
    // write_voxel_grid("bridge_test.stl", vertex, element);
    materials::HexDeformableBody<double> hex_def_body(linear_elasticity_material, vertex, dx, hex_mesh);

    // calculate constraints and external forces
    Eigen::Vector3d pmin = vertex.col(0), pmax = vertex.col(0);
    const int N = vertex.cols();
    for (int i = 0; i < N; ++i) {
        pmin = pmin.cwiseMin(vertex.col(i));
        pmax = pmax.cwiseMax(vertex.col(i));
    }
    std::vector<int> constraints;
    for (int i = 0; i < N; ++i)
        if (vertex.col(i)(0) == pmin(0) || vertex.col(i)(0) == pmax(0))
            constraints.push_back(i);

    Eigen::VectorXd F_ext(3 * N);
    F_ext.setZero();
    for (int i = 0; i < N; ++i)
        if (vertex.col(i)(2) == pmax(2))
            F_ext(i * 3 + 2) = -5000;

    compliance = hex_def_body.SolveCompliance(constraints, F_ext);
    num_voxels = element.cols();
}

int main(int argc, char *argv[])
{
    int N = 20000;
    std::vector<Eigen::Vector2d> p1_input; 
    std::vector<Eigen::VectorXd> p2_input;
    std::vector<Eigen::Vector2d> p3_input;
    std::vector<Eigen::Vector2d> p4_input;
    std::vector<Eigen::Vector3d> ec_input;

    p1_input.clear();
    p2_input.clear();
    p3_input.clear();
    ec_input.clear();

    // fix random seed to get repeatable results
    srand (1);

    for (int i = 0; i < N; i++) {
        p1_input.push_back(Eigen::Vector2d::Random()+Eigen::Vector2d::Ones());
        p2_input.push_back(Eigen::Vector4d::Random()+Eigen::Vector4d::Ones());
        p3_input.push_back(Eigen::Vector2d::Random()+Eigen::Vector2d::Ones());
        ec_input.push_back(Eigen::Vector3d::Random()+Eigen::Vector3d::Ones());
    }

    // test Q1: 2D convex hull
    std::vector<Eigen::Vector2d> p1_result = geometry::ConvexHull2D(p1_input);
    // print Q1 results
    std::ofstream file1;
    file1.open("q1_result.txt");
    file1 << "print Q1 test result" << std::endl;
    file1 << "Totol number of points: " << p1_result.size() << std::endl;
    for (int i = 0; i < p1_result.size(); i++) {
        file1 << "P" << i << ": " << std::endl;
        file1 << p1_result[i] << std::endl;
    }
    file1.close();

    // test Q2: naive Nd Pareto front
    std::vector<Eigen::VectorXd> p2_result = geometry::ParetoFrontNdNaive(p2_input);
    // print Q2 results
    std::ofstream file2;
    file2.open("q2_result.txt");
    file2 << "print Q2 test result" << std::endl;
    file2 << "Totol number of points: " << p2_result.size() << std::endl;
    for (int i = 0; i < p2_result.size(); i++) {
        file2 << "P" << i << ": " << std::endl;
        file2 << p2_result[i] << std::endl;
    }
    file2.close();

    // test Q3: fast 2d Pareto front
    std::vector<Eigen::Vector2d> p3_result = geometry::ParetoFront2D(p3_input);
    // print Q3 results
    std::ofstream file3;
    file3.open("q3_result.txt");
    file3 << "print Q3 test result" << std::endl;
    file3 << "Totol number of points: " << p3_result.size() << std::endl;
    for (int i = 0; i < p3_result.size(); i++) {
        file3 << "P" << i << ": " << std::endl;
        file3 << p3_result[i] << std::endl;
    }
    file3.close();

    // test EC: fast 3d pareto front
    std::vector<Eigen::Vector3d> ec_result = geometry::ParetoFront3D(ec_input);
    // print EC results
    std::ofstream fileEC1;
    fileEC1.open("ec_result.txt");
    fileEC1 << "print EC test result" << std::endl;
    fileEC1 << "Totol number of points: " << ec_result.size() << std::endl;
    for (int i = 0; i < ec_result.size(); i++) {
        fileEC1 << "P" << i << ": " << std::endl;
        fileEC1 << ec_result[i] << std::endl;
    }
    fileEC1.close();

    // Q4: implement the pipeline from design space to the performance space
    double compliance, num_voxels;

    // Once you debug the bridge example correct, comment the code before dash line
    // and uncomment the code after dash line to run the test on 121 bridges
    // SolvePerformance(PROJECT_SOURCE_DIR"/data/assignment6/bridge.stl", compliance, num_voxels);
    // SolvePerformance(PROJECT_SOURCE_DIR"/CSG/assn6_meshes/bridge_r_40_o_-25.stl", compliance, num_voxels);
    // std::cout << compliance << " " << num_voxels << std::endl;
    // -----------------------------------------------------------------------
    std::string base(PROJECT_SOURCE_DIR"/CSG/assn6_meshes/bridge");
    int radius_start = 30;
    int radius_end   = 40;
    int offset_start = -30;
    int offset_end   = -20;

    std::ofstream file4;
    file4.open("q4_result.txt");
    file4 << "print Q4 test result" << std::endl << std::endl;

    int count = 0;
    for (int r = radius_start; r <= radius_end; r++) {
        for (int o = offset_start; o <= offset_end; o++) {
            std::string mesh_name = base + "_r_" + std::to_string(r) + "_o_" + std::to_string(o) + ".stl";
            std::cout << "bridge_r_" + std::to_string(r) + "_o_" + std::to_string(o) << std::endl;
            SolvePerformance(mesh_name, compliance, num_voxels);
            p4_input.push_back(Eigen::Vector2d(compliance, num_voxels));
            file4 << "bridge_r_" + std::to_string(r) + "_o_" + std::to_string(o) + ".stl" << std::endl;
            file4 << "Compliance: " +  std::to_string(compliance) << std::endl;
            file4 << "Total mass: " + std::to_string(num_voxels) << std::endl << std::endl;
            count++;
        }
    }
    file4.close();

    std::vector<Eigen::Vector2d> p4_result = geometry::ParetoFront2D(p4_input);
    // print Q4 results
    file4.open("q4_pareto_front_result.txt");
    file4 << "print Q4 pareto front result" << std::endl;
    file4 << "Totol number of points: " << p4_result.size() << std::endl;
    for (int i = 0; i < p4_result.size(); i++) {
        file4 << "P" << i << ": " << std::endl;
        file4 << p4_result[i] << std::endl;
    }
    file4.close();

    // print data for visualization
    file4.open("result.txt");
    for (auto &v : p4_input)
        file4 << v(0) << " " << v(1) << std::endl;
    file4.close();

    return 0;
}
