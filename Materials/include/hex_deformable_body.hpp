// Tao Du
// taodu@csail.mit.edu
// Oct 12, 2016
#pragma once
#include "deformable_body.hpp"
#include "hexahedral_mesh.hpp"
#include "typedefs.hpp"
#include <cmath>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCholesky>
#include <string>
#include <sstream>
#include <fstream>

namespace materials {

    template<typename T>
    class HexDeformableBody : public DeformableBody<3, T> {
    public:
        HexDeformableBody(const Material<3, T>& material,
                          const Matrix3X<T>& initial_vertex_position,
                          const T density, const HexahedralMesh<T>& undeformed_hex_mesh)
                          : DeformableBody<3, T>(material, initial_vertex_position, undeformed_hex_mesh),
                            hex_size_((undeformed_hex_mesh.vertex(undeformed_hex_mesh.element(0)(0)) - undeformed_hex_mesh.vertex(undeformed_hex_mesh.element(0)(1))).norm()) {
                                Prepare();
                            }






        HexDeformableBody(const std::vector<std::reference_wrapper<const Material<3, T>>>& materials,
                          const std::vector<int>& material_id,
                          const Matrix3X<T>& initial_vertex_position,
                          const T density, const HexahedralMesh<T>& undeformed_hex_mesh)
                          : DeformableBody<3, T>(materials, material_id, initial_vertex_position, undeformed_hex_mesh),
                            hex_size_((undeformed_hex_mesh.vertex(undeformed_hex_mesh.element(0)(0)) - undeformed_hex_mesh.vertex(undeformed_hex_mesh.element(0)(1))).norm()) {
                                Prepare();
                            }

        HexDeformableBody(const HexDeformableBody& hex_body) : DeformableBody<3, T>(hex_body),
                                                               hex_size_(hex_body.hex_size_), K(hex_body.K) {}

        ~HexDeformableBody() {}


        typedef Eigen::Matrix<T, Eigen::Dynamic, 1> VectorXT;

        //TODO: Students should fill this out
	//vertices is a matrix of the current vertex positions (3 x n)        
        const Eigen::SparseMatrix<T> ComputeStiffnessMatrix(
                const Matrix3X<T>& vertices) const{
            std::vector<Eigen::Triplet<T>> triplet_list;
            const int vertex_num = static_cast<int>(this->vertex_position_.size() / 3);
           
            Eigen::SparseMatrix<T> K(vertex_num * 3, vertex_num * 3);

            // Beichen Li: the following code is based on uniform linear material
            const auto &elements = this->undeformed_mesh_.element();
            const int num_elements = elements.cols();

            // compute quadrature point coordinates (0-1)
            const T q[2] = {0.5 - 0.5 / std::sqrt(3.0), 0.5 + 0.5 / std::sqrt(3.0)};
            Eigen::Matrix<T, 3, 8> quad_pos;
            quad_pos << q[0], q[0], q[0], q[0], q[1], q[1], q[1], q[1],
                        q[0], q[0], q[1], q[1], q[0], q[0], q[1], q[1],
                        q[0], q[1], q[0], q[1], q[0], q[1], q[0], q[1];

            // compute dPdphi
            Eigen::Matrix<T, 3, 3> phi;
            const Material<3, T> &material = this->materials_[0];
            auto dPdphi = std::move(material.StressDifferential(phi));

            // compute local stiffness matrix
            const T inv_dx = 1.0 / hex_size_;
            const T coeff = hex_size_ * hex_size_ * hex_size_ / 8;
            Eigen::Matrix<T, 24, 24> stiffness_local;
            stiffness_local.setZero();
            for (int qi = 0; qi < 8; ++qi) {
                auto dphidx = std::move(DeformationGradientPartialx(quad_pos.col(qi), inv_dx));
                auto temp = std::move(dphidx.transpose() * dPdphi * dphidx * coeff);
                stiffness_local += (temp + temp.transpose()) * 0.5;
            }

            // reserve space for triplets
            triplet_list.reserve(24 * 24 * num_elements);

            // assemble global K matrix
            for (int ei = 0; ei < num_elements; ++ei)
                for (int i = 0; i < 8; ++i)
                    for (int j = 0; j < 8; ++j) {
                        int i_index = elements(i, ei) * 3;
                        int j_index = elements(j, ei) * 3;
                        for (int x = 0; x < 3; ++x)
                            for (int y = 0; y < 3; ++y)
                                triplet_list.emplace_back(Eigen::Triplet<T>(j_index + y, i_index + x,
                                    stiffness_local(j * 3 + y, i * 3 + x)));
                    }

            K.setFromTriplets(triplet_list.begin(), triplet_list.end());

            // Make sure K is symmetric.
            // K = (K + Eigen::SparseMatrix<T>(K.transpose())) / 2.0;
            return std::move(K);
        }

        // prepare for solving deformation
        void Prepare() {
            // calculate stiffness matrix
            this->K = std::move(ComputeStiffnessMatrix(this->undeformed_mesh_.vertex()));
        }

        // return stiffness matrix
        const Eigen::SparseMatrix<T> &StiffnessMatrix() { return this->K; }

        // test the difference between K and the reference answer
        T TestStiffnessMatrix(const std::string &file) {
            int num_vertices = this->undeformed_mesh_.NumOfVertex();

            // read reference K matrix from file
            std::ifstream fin(file.c_str());
            std::string line;
            std::vector<double> K_ref_in;
            K_ref_in.reserve(num_vertices * 3 * num_vertices * 3);
            while (std::getline(fin, line)) {
                std::stringstream sin(line);
                std::string valstr;
                while (std::getline(sin, valstr, ','))
                    K_ref_in.push_back(std::stod(valstr));
            }
            Eigen::Map<Eigen::MatrixXd> K_ref(K_ref_in.data(), num_vertices * 3, num_vertices * 3);

            // return maximum element difference
            Eigen::MatrixXd K_dense(this->K);
            return (K_ref - K_dense).lpNorm<Eigen::Infinity>();
        }

        // solve deformation given constraints and external forces
        // [params]
        //   vec: indices of fixed nodes
        //   F_ext: external forces
        VectorXT SolveDeformation(const std::vector<int> &vec, const VectorXT &F_ext) {
            const int num_constraints = vec.size();
            const int num_vertices = this->undeformed_mesh_.NumOfVertex();
            const int num_vertices_def = num_vertices - num_constraints;

            // generate permutation matrix
            std::vector<bool> con_flag(num_vertices, false);
            for (auto &v : vec)
                con_flag[v] = true;
            Eigen::VectorXi perm_indices(num_vertices * 3);
            int pos = 0;
            for (int i = 0; i < num_vertices; ++i)
                if (con_flag[i]) {
                    perm_indices.segment<3>(pos) = Eigen::Vector3i(i * 3, i * 3 + 1, i * 3 + 2);
                    pos += 3;
                }

            for (int i = 0; i < num_vertices; ++i)
                if (!con_flag[i]) {
                    perm_indices.segment<3>(pos) = Eigen::Vector3i(i * 3, i * 3 + 1, i * 3 + 2);
                    pos += 3;
                }
            Eigen::PermutationMatrix<Eigen::Dynamic> perm(perm_indices);

            // std::cout << perm.indices().transpose().eval() << std::endl;

            // remove rows and cols related to constrained nodes
            Eigen::SparseMatrix<T> K_perm(num_vertices * 3, num_vertices * 3);
            K_perm = this->K.twistedBy(perm.transpose());
            Eigen::SparseMatrix<T> K_con = std::move(K_perm.bottomRightCorner(num_vertices_def * 3, num_vertices_def * 3));

            // Eigen::MatrixXd K_con_dense(K_con);
            // std::cout << K_con_dense << std::endl;

            // solve for deformation
            Eigen::ConjugateGradient<Eigen::SparseMatrix<T>, Eigen::Upper|Eigen::Lower> cg;
            cg.compute(K_con);
            VectorXT U(num_vertices * 3);
            U.setZero();
            U.tail(num_vertices_def * 3) = cg.solve((perm.transpose() * F_ext).tail(num_vertices_def * 3));

            std::cout << "Iterations: " << cg.iterations() << ", Error: " << cg.error() << std::endl;

            return std::move(perm * U);
        }

        // solve deformation given constraints (fixed nodes indexed from 0 to num_constraints - 1) and external forces
        // [params]
        //   num_constraints: number of constrained nodes
        //   F_ext: external forces
        VectorXT SolveDeformation(int num_constraints, const VectorXT &F_ext) {
            const int num_vertices = this->undeformed_mesh_.NumOfVertex();
            const int num_vertices_def = num_vertices - num_constraints;

            Eigen::SparseMatrix<T> K_con = std::move(this->K.bottomRightCorner(num_vertices_def * 3, num_vertices_def * 3));

            // Eigen::MatrixXd K_con_dense(K_con);
            // std::cout << K_con_dense << std::endl;

            Eigen::ConjugateGradient<Eigen::SparseMatrix<T>, Eigen::Upper|Eigen::Lower> cg;
            cg.compute(K_con);
            VectorXT U(num_vertices * 3);
            U.setZero();
            U.tail(num_vertices_def * 3) = cg.solve(F_ext.tail(num_vertices_def * 3));

            std::cout << "Iterations: " << cg.iterations() << ", Error: " << cg.error() << std::endl;

            return std::move(U);
        }

        //return dphi (the deformation gradient) for a given voxel:
        //undeformed_vertex is a point in material space.
        //undeformed_cube are the vertices of an undeformed voxel
        //deformed_cube are the vertices of a the deformed voxel.
        static const Eigen::Matrix<T, 3, 3> DeformationGradient(
                const Vector3<T>& undeformed_vertex,
                const Eigen::Matrix<T, 3, 8>& undeformed_cube,
                const Eigen::Matrix<T, 3, 8>& deformed_cube){
            // Rename variables.
            const Vector3<T>& X = undeformed_vertex;
            const Eigen::Matrix<T, 3, 8>& X0 = undeformed_cube;
            const Eigen::Matrix<T, 3, 8>& x0 = deformed_cube;
            const T dx = X0(0, 4) - X0(0, 0);
            const T inv_dx = 1.0 / dx;
            const Vector3<T> offset = X - X0.col(0);
            const T rx = offset.x() / dx;
            const T ry = offset.y() / dx;
            const T rz = offset.z() / dx;
            Eigen::Matrix<T, 3, 3> F = Eigen::Matrix<T, 3, 3>::Zero();
            const T x_factor[2] = {1 - rx, rx};
            const T y_factor[2] = {1 - ry, ry};
            const T z_factor[2] = {1 - rz, rz};
            for (int i = 0; i < 2; ++i)
                for (int j = 0; j < 2; ++j)
                    for (int k = 0; k < 2; ++k) {
                        F.col(0) += x0.col(4 * i + 2 * j + k)
                                    * (i == 0 ? -inv_dx : inv_dx) * y_factor[j] * z_factor[k];
                        F.col(1) += x0.col(4 * i + 2 * j + k)
                                    * x_factor[i] * (j == 0 ? -inv_dx : inv_dx) * z_factor[k];
                        F.col(2) += x0.col(4 * i + 2 * j + k)
                                    * x_factor[i] * y_factor[j] * (k == 0 ? -inv_dx : inv_dx);
                    }
            return F;
        }




        //return dphi/dx for a given voxel:
        //undeformed_vertex is a point in material space.
        //undeformed_cube are the vertices of an undeformed voxel
        //deformed_cube are the vertices of a the deformed voxel.
        static const Eigen::Matrix<T, 9, 24> DeformationGradientPartialx(
                const Vector3<T>& undeformed_vertex,
                const Eigen::Matrix<T, 3, 8>& undeformed_cube,
                const Eigen::Matrix<T, 3, 8>& deformed_cube){

            const Vector3<T>& X = undeformed_vertex;
            const Eigen::Matrix<T, 3, 8>& X0 = undeformed_cube;
            Eigen::Matrix<T, 9, 24> Jacobian = MatrixX<T>::Zero(9, 24);
            const T dx = X0(0, 4) - X0(0, 0);
            const T inv_dx = 1.0 / dx;
            const Vector3<T> offset = X - X0.col(0);
            const T rx = offset.x() / dx;
            const T ry = offset.y() / dx;
            const T rz = offset.z() / dx;
            Eigen::Matrix<T, 3, 3> F = Eigen::Matrix<T, 3, 3>::Zero();
            const T x_factor[2] = { 1 - rx, rx };
            const T y_factor[2] = { 1 - ry, ry };
            const T z_factor[2] = { 1 - rz, rz };
            for (int i = 0; i < 2; ++i)
                for (int j = 0; j < 2; ++j)
                    for (int k = 0; k < 2; ++k) {
                        const int index = 4 * i + 2 * j + k;
                        const T scale_first_column = (i == 0 ? -inv_dx : inv_dx)
                                                          * y_factor[j] * z_factor[k];
                        Jacobian(0, 3 * index) += scale_first_column;
                        Jacobian(1, 3 * index + 1) += scale_first_column;
                        Jacobian(2, 3 * index + 2) += scale_first_column;
                        const T scale_second_column = x_factor[i]
                                                           * (j == 0 ? -inv_dx : inv_dx) * z_factor[k];
                        Jacobian(3, 3 * index) += scale_second_column;
                        Jacobian(4, 3 * index + 1) += scale_second_column;
                        Jacobian(5, 3 * index + 2) += scale_second_column;
                        const T scale_third_column = x_factor[i] * y_factor[j]
                                                          * (k == 0 ? -inv_dx : inv_dx);
                        Jacobian(6, 3 * index) += scale_third_column;
                        Jacobian(7, 3 * index + 1) += scale_third_column;
                        Jacobian(8, 3 * index + 2) += scale_third_column;
                    }
            return Jacobian;



        }

        // Beichen Li: DeformationGradientPartialx for a quadrature point (represented in coefficient matrix) in uniform grid
        // [params]
        //   r: relative coordinates of quadrature point (in [0, 1]^3)
        //   inv_dx: 1.0 / size of each hexahedron element
        static const Eigen::Matrix<T, 9, 24> DeformationGradientPartialx(const Vector3<T> &r, const T inv_dx){
            // Rename variables.
            const T rx = r(0);
            const T ry = r(1);
            const T rz = r(2);
            const T x_factor[2] = {1 - rx, rx};
            const T y_factor[2] = {1 - ry, ry};
            const T z_factor[2] = {1 - rz, rz};

            // compute coefficient matrix
            Eigen::Matrix<T, 3, 8> coeff;
            int row_index;
            for (int i = 0; i < 2; ++i)
                for (int j = 0; j < 2; ++j)
                    for (int k = 0; k < 2; ++k) {
                        int index = i * 4 + j * 2 + k;
                        coeff(0, index) = (i == 0 ? -inv_dx : inv_dx) * y_factor[j] * z_factor[k];
                        coeff(1, index) = x_factor[i] * (j == 0 ? -inv_dx : inv_dx) * z_factor[k];
                        coeff(2, index) = x_factor[i] * y_factor[j] * (k == 0 ? -inv_dx : inv_dx);
                    }

            // compute Jacobian
            Eigen::Matrix<T, 9, 24> Jacobian;
            Jacobian.setZero();
            for (int i = 0; i < 8; ++i)
                for (int j = 0; j < 3; ++j)
                    for (int k = 0; k < 3; ++k)
                        Jacobian(j * 3 + k, i * 3 + k) = coeff(j, i);

            return std::move(Jacobian);
        }

        virtual void getInitialNodes(Matrix3X<T>& initial_nodes){
            initial_nodes = this->undeformed_mesh_.vertex();
        }

        void GaussIntegrationFactorCheck() {
            Eigen::Matrix<T, 3, 8> coords;
            coords << 0, 0, 0, 0, 1, 1, 1, 1,
                      0, 0, 1, 1, 0, 0, 1, 1,
                      0, 1, 0, 1, 0, 1, 0, 1;
            coords = coords * GaussIntegrationFactor();

            Eigen::Matrix<T, 3, 8> quad_pos;
            const T q[2] = {0.5 - 0.5 / std::sqrt(3.0), 0.5 + 0.5 / std::sqrt(3.0)};
            quad_pos << q[0], q[0], q[0], q[0], q[1], q[1], q[1], q[1],
                        q[0], q[0], q[1], q[1], q[0], q[0], q[1], q[1],
                        q[0], q[1], q[0], q[1], q[0], q[1], q[0], q[1];

            std::cout << coords << std::endl;
            std::cout << quad_pos << std::endl;
        }

    private:
        HexDeformableBody& operator=(const HexDeformableBody&);


        //TODO: Studnets fill this function out
        static const Eigen::Matrix<T, 8, 8> GaussIntegrationFactor() {
            // \int_{-1}^{1} f(x) dx \approx f(-1/sqrt(3)) + f(1/sqrt(3)).
            Eigen::Matrix<T, 8, 8> X0_coeff = Eigen::MatrixXd::Zero(8, 8);
            
            const T e_pos = 1.0 + 1.0 / std::sqrt(3.0);
            const T e_neg = 1.0 - 1.0 / std::sqrt(3.0);
            const T e[8] = {
                e_neg * e_neg * e_neg / 8,
                e_neg * e_neg * e_pos / 8,
                e_neg * e_pos * e_neg / 8,
                e_neg * e_pos * e_pos / 8,
                e_pos * e_neg * e_neg / 8,
                e_pos * e_neg * e_pos / 8,
                e_pos * e_pos * e_neg / 8,
                e_pos * e_pos * e_pos / 8
            };

            for (int j = 0; j < 8; ++j)
                for (int i = 0; i < 8; ++i)
                    X0_coeff(i, j) = e[i ^ j ^ 7];

            return X0_coeff;
        }


        const T hex_size_;

        // Beichen Li: ready to solve for deformation
        Eigen::SparseMatrix<T> K;
    };

}
