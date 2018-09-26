#pragma once

#include "read_stl.hpp"
#include "BasicGeometry.hpp"
#include <math.h>
#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
#include <Eigen/Core>
#include <unordered_set>
#include <fstream>
#include <ctime>

namespace mesh {

    template<typename T>
    class Voxelizer {
    public:
        Voxelizer(const std::string& stl_file_name, const T dx)
            : _dx(dx) {
            // Randomness.
            srand(static_cast<unsigned>(time(0)));
            // Load triangles from the stl file.
            std::vector<Vector3<T>> normals;
            if (!ReadSTL(stl_file_name, _triangles, normals)) {
                std::cout << "ERROR: cannot read " << stl_file_name << std::endl;
                return;
            }
            // Compute the bounding box of _triangle and save the results into _pmin.
            _pmin = _triangles[0][0];
            Vector3<T> pmax = _triangles[0][0];
            for (const auto& triangle : _triangles)
                for (const auto& v : triangle) {
                    _pmin = _pmin.cwiseMin(v);
                    pmax = pmax.cwiseMax(v);
                }
            for (int i = 0; i < 3; ++i) {
                _pmin[i] -= _dx;
                pmax[i] += _dx;
            }
            // Compute the number of voxels along each direction.
            for (int i = 0; i < 3; ++i)
                _nvoxel[i] = static_cast<int>((pmax[i] - _pmin[i]) / _dx) + 1;
            // Initialize the voxel array.
            _voxels = std::vector<std::vector<std::vector<bool>>>(_nvoxel.x(),
                std::vector<std::vector<bool>>(_nvoxel.y(),
                    std::vector<bool>(_nvoxel.z(), false)));
        }

        const Vector3<T> pmin() const { return _pmin; }
        const T dx() const { return _dx; }
        const Vector3<T> pmax() const { return _pmin + Vector3<T>(_nvoxel.x(), _nvoxel.y(), _nvoxel.z()) * _dx; }
        const Vector3<int> voxel_num() const { return _nvoxel; }

        void BasicVoxelization() {
            /* Assignment 2, Part 2.1. */
            /* Implement your code here. */
            // Fill the _voxels array with the correct flag.
            const int nx = _nvoxel[0], ny = _nvoxel[1], nz = _nvoxel[2];
            for (int i = 0; i < nx; ++i)
                for (int j = 0; j < ny; ++j)
                    for (int k = 0; k < nz; ++k)
                        _voxels[i][j][k] = false;

            // Store all intersection points (in parameter t)
            std::vector<T> vec;

            // Shoot an upright ray from each grid cell on the bottom plane
            T curx = _pmin(0) + 0.5 * _dx;
            T curz = _pmin(2);
            Vector3<T> dir(0, 0, 1);
            for (int i = 0; i < nx; ++i, curx += _dx) {
                T cury = _pmin(1) + 0.5 * _dx;
                for (int j = 0; j < ny; ++j, cury += _dx) {
                    // Ray-triangle intersection
                    Vector3<T> origin(curx, cury, curz);
                    vec.clear();
                    for (const auto &_tri : _triangles) {
                        geometry::Triangle<T> tri(_tri[0], _tri[1], _tri[2]);
                        T t = tri.IntersectRay(origin, dir);
                        if (t >= 0.0)
                            vec.push_back(t);
                    }
                    if (vec.empty()) continue;

                    // Sort intersection points and eliminate duplication
                    Deduplicate(vec);
                    // printf("(%d, %d) - pts: %d\n", i, j, (int)vec.size());

                    // Fill _voxels array
                    for (int idx = 0; idx + 1 < (int)vec.size(); idx += 2) {
                        int start = std::max((int)std::ceil(vec[idx] / _dx - 0.5 - 1e-6), 0);
                        int end = std::min((int)std::floor(vec[idx + 1] / _dx - 0.5 + 1e-6), nz - 1);
                        for (int k = start; k <= end; ++k)
                            _voxels[i][j][k] = true;
                    }
                }
            }
        }

        void AdvancedVoxelization() {
            /* Assignment 2, Part 2.2. */
            /* Implement your code here. */
            // Fill the _voxels array with the correct flag.
            const int nx = _nvoxel[0], ny = _nvoxel[1], nz = _nvoxel[2];
            for (int i = 0; i < nx; ++i)
                for (int j = 0; j < ny; ++j)
                    for (int k = 0; k < nz; ++k)
                        _voxels[i][j][k] = false;

            // Store intersection points for all rays
            std::vector<T> vec[nx][ny];


        }

        void AdvancedVoxelizationWithApproximation() {
            /* Assignment 2, Part 2.3. */
            /* Implement your code here. */
            // Fill the _voxels array with the correct flag.
            const int nx = _nvoxel[0], ny = _nvoxel[1], nz = _nvoxel[2];
            for (int i = 0; i < nx; ++i)
                for (int j = 0; j < ny; ++j)
                    for (int k = 0; k < nz; ++k)
                        _voxels[i][j][k] = false;
        }

        void WriteVoxelToMesh(const std::string& stl_file_name) const {
            const int nx = _nvoxel[0], ny = _nvoxel[1], nz = _nvoxel[2];
            std::vector<std::vector<Vector3<int>>> faces;
            std::vector<Vector3<int>> corners({
                Vector3<int>(0, 0, 0),
                Vector3<int>(0, 0, 1),
                Vector3<int>(0, 1, 0),
                Vector3<int>(0, 1, 1),
                Vector3<int>(1, 0, 0),
                Vector3<int>(1, 0, 1),
                Vector3<int>(1, 1, 0),
                Vector3<int>(1, 1, 1)
            });
            for (int i = 0; i < nx; ++i)
                for (int j = 0; j < ny; ++j)
                    for (int k = 0; k < nz; ++k) {
                        if (!_voxels[i][j][k]) continue;
                        // Check -x direction.
                        Vector3<int> cmin(i, j, k);
                        if (i == 0 || !_voxels[i - 1][j][k]) {
                            faces.push_back({ cmin + corners[0], cmin + corners[1], cmin + corners[3] });
                            faces.push_back({ cmin + corners[0], cmin + corners[3], cmin + corners[2] });
                        }
                        if (i == nx - 1 || !_voxels[i + 1][j][k]) {
                            faces.push_back({ cmin + corners[4], cmin + corners[6], cmin + corners[7] });
                            faces.push_back({ cmin + corners[4], cmin + corners[7], cmin + corners[5] });
                        }
                        if (j == 0 || !_voxels[i][j - 1][k]) {
                            faces.push_back({ cmin + corners[0], cmin + corners[4], cmin + corners[5] });
                            faces.push_back({ cmin + corners[0], cmin + corners[5], cmin + corners[1] });
                        }
                        if (j == ny - 1 || !_voxels[i][j + 1][k]) {
                            faces.push_back({ cmin + corners[2], cmin + corners[3], cmin + corners[7] });
                            faces.push_back({ cmin + corners[2], cmin + corners[7], cmin + corners[6] });
                        }
                        if (k == 0 || !_voxels[i][j][k - 1]) {
                            faces.push_back({ cmin + corners[0], cmin + corners[2], cmin + corners[6] });
                            faces.push_back({ cmin + corners[0], cmin + corners[6], cmin + corners[4] });
                        }
                        if (k == nz - 1 || !_voxels[i][j][k + 1]) {
                            faces.push_back({ cmin + corners[5], cmin + corners[7], cmin + corners[3] });
                            faces.push_back({ cmin + corners[5], cmin + corners[3], cmin + corners[1] });
                        }
                    }
            std::ofstream fout(stl_file_name);
            fout << "solid vcg" << std::endl;
            for (const auto& f : faces) {
                std::vector<Vector3<T>> p;
                for (const auto& fi : f) {
                    Vector3<T> v = _pmin + fi.cast<T>() * _dx;
                    p.push_back(v);
                }
                const Vector3<T> n = (p[1] - p[0]).cross(p[2] - p[1]).normalized();
                fout << "  facet normal " << n.x() << " " << n.y() << " " << n.z() << std::endl;
                fout << "    outer loop" << std::endl;
                for (const auto& v : p) {
                    fout << "      vertex " << v.x() << " " << v.y() << " " << v.z() << std::endl;
                }
                fout << "    endloop" << std::endl;
                fout << "  endfacet" << std::endl;
            }
            fout << "endsolid vcg" << std::endl;
        }

    private:
        std::vector<std::vector<Vector3<T>>> _triangles;
        T _dx;  // The size of each voxel.
        Vector3<T> _pmin;    // The min and max corner of the bounding box.
        Eigen::Vector3i _nvoxel;   // The number of voxels along each direction.
        std::vector<std::vector<std::vector<bool>>> _voxels;   // True <-> voxel is occupied.

        // Deduplication of floating-point vector (with error tolerance)
        void Deduplicate(std::vector<T> &vec) {
            std::sort(vec.begin(), vec.end());
            int pos = 0;
            for (int i = 1; i < (int)vec.size(); ++i)
                if (vec[i] - vec[pos] > 1e-6) // Ascending order as default
                    vec[++pos] = vec[i];
            vec.erase(vec.begin() + pos + 1, vec.end());
        }
    };

}
