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




#include <stdio.h>
#include <iostream>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCholesky>
#include <fstream>
#include <set>

#include "voxelizer.hpp"

#define STR(x)  #x


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


int main(int argc, char *argv[])
{

    const int dim = 3;

    //youngs =  10^7, poisson ratio = 0.45 - similar to that of silicone rubber
    materials::LinearElasticityMaterial<dim, double> linear_elasticity_material(10000000, 0.45);

    const size_t num_x_vertices = 5; //number of nodes (vertices) in x dimension of our mesh
    const size_t num_y_vertices = 3; //number of nodes (vertices) in y dimension of our mesh
    const size_t num_z_vertices = 2; //number of nodes (vertices) in z dimension of our mesh

    const size_t num_vertices = num_x_vertices * num_y_vertices * num_z_vertices;

    const double spacing = 0.05; //size of each edge, 5 cm
    materials::HexahedralMesh<double> hex_mesh =
            materials::HexahedralMeshCuboid<double>(Eigen::Vector3d::Zero(),
                    Eigen::Vector3i(num_x_vertices, num_y_vertices, num_z_vertices), spacing);

    //What does our undeformed mesh look like?
    write_voxel_grid("test.stl", hex_mesh.vertex(), hex_mesh.element());



    materials::HexDeformableBody<double> hex_def_body(linear_elasticity_material, hex_mesh.vertex(), 0.4, hex_mesh);

    //TODO: Students take it from here!

    /// Simulation on uniform grid

    // test stiffness matrix
    // std::cout << "Maximum error: " <<
    //     hex_def_body.TestStiffnessMatrix("../ComputationalFabrication/data/assignment4/test.txt") << std::endl;
    
    // external forces
    const int nx = num_x_vertices;
    const int ny = num_y_vertices;
    const int nz = num_z_vertices;
    Eigen::VectorXd F_ext(3 * num_vertices);
    F_ext.setZero();
    for (int y = 0; y < ny; ++y)
        F_ext(((nx - 1) * ny + y) * nz * 3 + 2) = -100;

    // solve deformation
    Eigen::VectorXd U = std::move(hex_def_body.SolveDeformation(ny * nz, F_ext));
    // std::cout << U.transpose().eval() << std::endl;

    // show deformed mesh
    Eigen::Map<Eigen::MatrixXd> U_mat(U.data(), 3, num_vertices);
    write_voxel_grid("test_deformed.stl", hex_mesh.vertex() + U_mat, hex_mesh.element());

    /// Simulation on customized mesh

    // import mesh from file
    mesh::Voxelizer<double> voxelizer("../ComputationalFabrication/data/assignment2/fandisk.stl", 0.4);
    voxelizer.AdvancedVoxelization();
    materials::HexahedralMesh<double> hex_mesh_c = std::move(voxelizer.ConvertToHexMesh());
    
    // write initial mesh
    write_voxel_grid("test_c.stl", hex_mesh_c.vertex(), hex_mesh_c.element());
    // std::cout << "Vertices: " << hex_mesh_c.vertex().cols() << std::endl;
    // std::cout << "Elements: " << hex_mesh_c.element().cols() << std::endl;
    // std::cout << hex_mesh_c.element() << std::endl;

    materials::HexDeformableBody<double> hex_def_body_c(linear_elasticity_material, hex_mesh_c.vertex(), 0.4, hex_mesh_c);

    // Eigen::SparseMatrix<double> K_c = hex_def_body_c.StiffnessMatrix();
    // Eigen::MatrixXd K_dense_c(K_c);
    // std::cout << K_dense_c << std::endl;

    // add constraints (set points with minimal z coordinates to be fixed)
    const auto &vertex_c = hex_mesh_c.vertex();
    const int num_vertices_c = vertex_c.cols();
    Eigen::Vector3d pmin = vertex_c.col(0);
    Eigen::Vector3d pmax = vertex_c.col(0);
    for (int i = 1; i < num_vertices_c; ++i) {
        pmin = std::move(pmin.cwiseMin(vertex_c.col(i)));
        pmax = std::move(pmax.cwiseMax(vertex_c.col(i)));
    }
    std::vector<int> vertices_con;
    vertices_con.reserve(32);
    for (int i = 0; i < num_vertices_c; ++i)
        if (vertex_c(2, i) == pmin(2))
            vertices_con.push_back(i);

    std::cout << "Constraints: " << vertices_con.size() << std::endl;

    // external forces (apply force in z direction to points with maximal x coordinates)
    Eigen::VectorXd F_ext_c(3 * num_vertices_c);
    F_ext_c.setZero();
    for (int i = 0; i < num_vertices_c; ++i)
        if (vertex_c(1, i) == pmax(1))
            F_ext_c(i * 3 + 2) = 5000;

    // solve deformation
    Eigen::VectorXd U_c = std::move(hex_def_body_c.SolveDeformation(vertices_con, F_ext_c));
    std::cout << "Maximum displacement: " << U_c.lpNorm<Eigen::Infinity>() << std::endl;

    // show deformed mesh
    Eigen::Map<Eigen::MatrixXd> U_mat_c(U_c.data(), 3, num_vertices_c);
    write_voxel_grid("test_c_deformed.stl", vertex_c + U_mat_c, hex_mesh_c.element());


    std::cout << "Done with OpenFab!  Have a squishy day." << std::endl;

}
