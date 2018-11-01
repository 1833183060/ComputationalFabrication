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

    Eigen::SparseMatrix<double> K = std::move(hex_def_body.ComputeStiffnessMatrix(hex_mesh.vertex()));
    
    // Beichen Li: test K matrix
    // Eigen::Matrix<double, num_vertices * 3, num_vertices * 3> K_dense(K);
    // std::ofstream file("test.txt", std::ios::out);
    // file << K_dense << std::endl;
    // file.close();

    const int nx = num_x_vertices;
    const int ny = num_y_vertices;
    const int nz = num_z_vertices;
    
    // external forces
    Eigen::VectorXd F_ext(3 * num_vertices);
    F_ext.setZero();
    for (int y = 0; y < ny; ++y)
        F_ext(((nx - 1) * ny + y) * nz * 3 + 2) = -100;

    // add constraints
    const int num_vertices_con = num_vertices - ny * nz;
    Eigen::SparseMatrix<double> K_c = std::move(K.bottomRightCorner(3 * num_vertices_con, 3 * num_vertices_con));
    Eigen::VectorXd F_ext_c = std::move(F_ext.tail(3 * num_vertices_con));
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Upper|Eigen::Lower> cg;
    cg.compute(K_c);

    // solve deformation
    Eigen::VectorXd U(3 * num_vertices);
    U.setZero();
    U.tail(3 * num_vertices_con) = std::move(cg.solve(F_ext_c));
    std::cout << "CG iteration: " << cg.iterations() << ", error = " << cg.error() << std::endl;
    std::cout << U.transpose().eval() << std::endl;

    // show deformed mesh
    Eigen::Map<Eigen::MatrixXd> U_mat(U.data(), 3, num_vertices);
    write_voxel_grid("test_deformed.stl", hex_mesh.vertex() + U_mat, hex_mesh.element());

    std::cout << "Done with OpenFab!  Have a squishy day." << std::endl;

}
