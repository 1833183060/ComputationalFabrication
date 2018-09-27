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

            // Allocate triangle list
            nTriangles = _triangles.size();
            _tris = (geometry::Triangle<T> *)std::malloc(sizeof(geometry::Triangle<T>) * nTriangles);
            for (int i = 0; i < nTriangles; ++i)
                new (&_tris[i]) geometry::Triangle<T>(_triangles[i][0], _triangles[i][1], _triangles[i][2]);
        }

        ~Voxelizer() {
            std::free(_tris);
            _tris = NULL;
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
            std::vector<T> rays[nx][ny];

            // Shoot an upright ray from each grid cell on the bottom plane
            T curx = _pmin(0) + 0.5 * _dx;
            T curz = _pmin(2);
            Vector3<T> dir(0, 0, 1);
            for (int i = 0; i < nx; ++i, curx += _dx) {
                T cury = _pmin(1) + 0.5 * _dx;
                for (int j = 0; j < ny; ++j, cury += _dx) {
                    std::vector<T> &vec = rays[i][j];

                    // Ray-triangle intersection
                    Vector3<T> origin(curx, cury, curz);
                    for (int it = 0; it < nTriangles; ++it) {
                        T t = _tris[it].IntersectRay(origin, dir);
                        if (t >= 0.0)
                            vec.push_back(t);
                    }
                    if (vec.empty()) continue;

                    // Sort intersection points and eliminate duplication
                    int vecSize = Deduplicate(vec);

                    // Fill _voxels array
                    for (int idx = 0; idx + 1 < vecSize; idx += 2) {
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
            std::vector<T> rays[nx][ny];

            // Check each triangle for grid cell coverage
            for (int it = 0; it < nTriangles; ++it) {
                // Calculate bounding box for each triangle
                Vector3<T> tmin = std::move(_tris[it].bmin());
                Vector3<T> tmax = std::move(_tris[it].bmax());
                if (tmax(0) - tmin(0) < 1e-6 || tmax(1) - tmin(1) < 1e-6)
                    continue;
                int xs = std::max((int)std::ceil((tmin(0) - _pmin(0)) / _dx - 0.5 - 1e-6), 0);
                int xe = std::min((int)std::floor((tmax(0) - _pmin(0)) / _dx - 0.5 + 1e-6), nx - 1);
                int ys = std::max((int)std::ceil((tmin(1) - _pmin(1)) / _dx - 0.5 - 1e-6), 0);
                int ye = std::min((int)std::floor((tmax(1) - _pmin(1)) / _dx - 0.5 + 1e-6), ny - 1);

                // Shoot rays and check intersection
                T curx, cury;
                curx = _pmin(0) + (xs + 0.5) * _dx;
                for (int i = xs; i <= xe; ++i, curx += _dx) {
                    cury = _pmin(1) + (ys + 0.5) * _dx;
                    for (int j = ys; j <= ye; ++j, cury += _dx) {
                        T t = _tris[it].IntersectRay(Vector3<T>(curx, cury, _pmin(2)), Vector3<T>(0, 0, 1));
                        if (t >= 0.0)
                            rays[i][j].push_back(t);
                    }
                }
            }

            // Fill _voxels array
            for (int i = 0; i < nx; ++i)
                for (int j = 0; j < ny; ++j) {
                    std::vector<T> &vec = rays[i][j];
                    int vecSize = Deduplicate(vec);
                    for (int idx = 0; idx + 1 < vecSize; idx += 2) {
                        int start = std::max((int)std::ceil(vec[idx] / _dx - 0.5 - 1e-6), 0);
                        int end = std::min((int)std::floor(vec[idx + 1] / _dx - 0.5 + 1e-6), nz - 1);
                        for (int k = start; k <= end; ++k)
                            _voxels[i][j][k] = true;
                    }
                }
        }

        void AdvancedVoxelizationWithApproximation() {
            /* Assignment 2, Part 2.3. */
            /* Implement your code here. */
            // Initialization
            const int nx = _nvoxel[0], ny = _nvoxel[1], nz = _nvoxel[2];
            const Vector3<T> _pmin_backup = _pmin;
            auto &&_voxels_backup = std::move(_voxels);

            // Statistics of occupation
            char stat[nx][ny][nz] = {0};

            // History of rotation
            Eigen::Quaternion<T> q0 = Eigen::Quaternion<T>::Identity();

            // Main loop: rotate mesh by arbitrary angle and recalculate occupancy information
            // New occupancy values will be integrated into statistics with nearest neighbor sampling.
            int iter = 11;
            while (iter--) {
                // Mesh rotation
                Eigen::Quaternion<T> q = std::move(Eigen::Quaternion<T>::UnitRandom());
                for (int i = 0; i < nTriangles; ++i)
                    _tris[i].rotate(q);

                // Record rotation history
                q0 = q * q0;

                // Calculate bounding box
                _pmin = _tris[0].vertices(0);
                Vector3<T> pmax = _pmin;
                for (int i = 0; i < nTriangles; ++i)
                    for (int j = 0; j < 3; ++j) {
                        _pmin = _pmin.cwiseMin(_tris[i].vertices(j));
                        pmax = pmax.cwiseMax(_tris[i].vertices(j));
                    }
                _pmin -= Vector3<T>(_dx, _dx, _dx);
                pmax += Vector3<T>(_dx, _dx, _dx);

                // Allocate occupancy grid
                for (int i = 0; i < 3; ++i)
                    _nvoxel[i] = (int)std::ceil((pmax(i) - _pmin(i)) / _dx + 1e-6);
                _voxels = std::move(std::vector<std::vector<std::vector<bool>>>(_nvoxel[0],
                    std::vector<std::vector<bool>>(_nvoxel[1],
                        std::vector<bool>(_nvoxel[2], false))));

                // Voxelization
                AdvancedVoxelization();

                // Collect result
                for (int i = 0; i < _nvoxel[0]; ++i)
                    for (int j = 0; j < _nvoxel[1]; ++j)
                        for (int k = 0; k < _nvoxel[2]; ++k) {
                            // Convert coordinates back to original occupancy grid
                            Vector3<T> coord = std::move(Vector3<T>(i + 0.5, j + 0.5, k + 0.5) * _dx);
                            Vector3<T> p = std::move(q0.toRotationMatrix().transpose() * (_pmin + coord));

                            int px = (int)std::round((p(0) - _pmin_backup(0)) / _dx - 0.5);
                            int py = (int)std::round((p(1) - _pmin_backup(1)) / _dx - 0.5);
                            int pz = (int)std::round((p(2) - _pmin_backup(2)) / _dx - 0.5);

                            if (px >= 0 && px < nx && py >= 0 && py < ny && pz >= 0 && pz < nz)
                                stat[px][py][pz] += _voxels[i][j][k] ? 1 : -1;
                        }
            }

            // Compute final representation
            _nvoxel[0] = nx, _nvoxel[1] = ny, _nvoxel[2] = nz;
            _voxels = _voxels_backup;
            for (int i = 0; i < nx; ++i)
                for (int j = 0; j < ny; ++j)
                    for (int k = 0; k < nz; ++k)
                        _voxels[i][j][k] = stat[i][j][k] > 0;

            // Restore original values
            _pmin = _pmin_backup;
            for (int i = 0; i < nTriangles; ++i)
                new (&_tris[i]) geometry::Triangle<T>(_triangles[i][0], _triangles[i][1], _triangles[i][2]);
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

        // Triangle object list
        int nTriangles;
        geometry::Triangle<T> *_tris = NULL;

        // Deduplication of floating-point vector (with error tolerance)
        int Deduplicate(std::vector<T> &vec) {
            std::sort(vec.begin(), vec.end());
            int pos = 0;
            for (int i = 1; i < (int)vec.size(); ++i)
                if (vec[i] - vec[pos] > 1e-6) // Ascending order as default
                    vec[++pos] = vec[i];
            // vec.erase(vec.begin() + pos + 1, vec.end());
            return pos + 1;
        }
    };

}
