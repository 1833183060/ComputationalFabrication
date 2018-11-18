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

// get performance metrics
std::pair<double, int> getPerformance(const std::string &filename,
    const materials::Material<3, double> &material, const double dx) {

    mesh::Voxelizer<double> voxelizer(filename, dx);
    voxelizer.AdvancedVoxelization();
    materials::HexahedralMesh<double> hex_mesh = std::move(voxelizer.ConvertToHexMesh());
    const auto &vertex = hex_mesh.vertex();
    const auto &element = hex_mesh.element();
    materials::HexDeformableBody<double> hex_def_body(material, vertex, dx, hex_mesh);

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

    // solve linear equations for deformation
    // auto &&U = hex_def_body.SolveDeformation(constraints, F_ext);
    // Eigen::Map<Eigen::MatrixXd> U_mat(U.data(), 3, N);
    // write_voxel_grid("bridge_deform.stl", vertex + U_mat, element);

    return std::make_pair(hex_def_body.SolveCompliance(constraints, F_ext), element.cols());
}

int main(int argc, char *argv[])
{
    int N = 20000;
    std::vector<Eigen::Vector2d> p1_input; 
    std::vector<Eigen::VectorXd> p2_input;
    std::vector<Eigen::Vector2d> p3_input;
    std::vector<Eigen::Vector3d> ec1_input;
    std::vector<Eigen::Vector3d> ec2_input;

    p1_input.clear();
    p2_input.clear();
    p3_input.clear();
    ec1_input.clear();
    ec2_input.clear();

    // fix random seed to get repeatable results
    srand (1);

    for (int i = 0; i < N; i++) {
        p1_input.push_back(Eigen::Vector2d::Random()+Eigen::Vector2d::Ones());
        p2_input.push_back(Eigen::Vector4d::Random()+Eigen::Vector4d::Ones());
        p3_input.push_back(Eigen::Vector2d::Random()+Eigen::Vector2d::Ones());
        ec1_input.push_back(Eigen::Vector3d::Random()+Eigen::Vector3d::Ones());
        ec2_input.push_back(Eigen::Vector3d::Random()+Eigen::Vector3d::Ones());
    }

    // test Q1: 2D convex hull
    std::vector<Eigen::Vector2d> p1_result = geometry::ConvexHull2D(p1_input);
    // print Q1 results
    std::ofstream file1;
    file1.open("q1_result2.txt");
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
    file2.open("q2_result2.txt");
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
    file3.open("q3_result2.txt");
    file3 << "print Q3 test result" << std::endl;
    file3 << "Totol number of points: " << p3_result.size() << std::endl;
    for (int i = 0; i < p3_result.size(); i++) {
        file3 << "P" << i << ": " << std::endl;
        file3 << p3_result[i] << std::endl;
    }
    file3.close();

    // test EC1: fast 3d convex hull
    std::vector<Eigen::Vector3d> ec1_result = geometry::ConvexHull3D(ec1_input);
    // print EC1 results
    // std::ofstream fileEC1;
    // fileEC1.open("ec1_result.txt");
    // fileEC1 << "print EC1 test result" << std::endl;
    // fileEC1 << "Totol number of points: " << ec1_result.size() << std::endl;
    // for (int i = 0; i < ec1_result.size(); i++) {
    //     fileEC1 << "P" << i << ": " << std::endl;
    //     fileEC1 << ec1_result[i] << std::endl;
    // }
    // fileEC1.close();

    // test EC2: fast 3d pareto front
    std::vector<Eigen::Vector3d> ec2_result = geometry::ParetoFront3D(ec2_input);
    // print EC2 results
    std::ofstream fileEC2;
    fileEC2.open("ec2_result.txt");
    fileEC2 << "print EC2 test result" << std::endl;
    fileEC2 << "Totol number of points: " << ec2_result.size() << std::endl;
    for (int i = 0; i < ec2_result.size(); i++) {
        fileEC2 << "P" << i << ": " << std::endl;
        fileEC2 << ec2_result[i] << std::endl;
    }
    fileEC2.close();

    // Q4: implement the pipeline from design space to the performance space
    const int dim = 3;
    materials::LinearElasticityMaterial<dim, double> linear_elasticity_material(10000000, 0.45);
    // 
    // voxelize bridge
    std::string stl_name(PROJECT_SOURCE_DIR"/CSG/assn6_meshes/bridge.stl");
    double dx = 0.25;

    // test sample bridge
    // auto &&perf_test = getPerformance(stl_name, linear_elasticity_material, dx);
    // std::cout << perf_test.first << " " << perf_test.second << std::endl;
    // return 0;

    std::string dir(PROJECT_SOURCE_DIR"/CSG/assn6_meshes/");
    std::vector<Eigen::Vector2d> gamut;
    for (int r = 30; r <= 40; ++r)
        for (int o = -30; o <= -20; ++o) {
            std::stringstream ss;
            ss << dir << "bridge_r_" << r << "_o_" << o << ".stl";
            std::cout << "r = " << r / 10.0 << ", o = " << o / 10.0 << std::endl;

            auto &&perf = getPerformance(ss.str(), linear_elasticity_material, dx);
            gamut.push_back(Eigen::Vector2d(perf.first, perf.second));
        }

    auto &&pareto_front = geometry::ParetoFront2D(gamut);
    std::cout << "Collected " << gamut.size() << " samples" << std::endl;
    
    // save gamut
    std::ofstream fout("gamut.txt");
    for (auto &p : gamut)
        fout << p(0) << " " << p(1) << std::endl;
    fout.close();

    // Pareto front
    std::cout << "Pareto front: " << pareto_front.size() << " samples" << std::endl;
    for (auto &p : pareto_front)
        std::cout << p(0) << " " << p(1) << std::endl;

    return 0;
}
