#pragma once
#include "tri_mesh.hpp"
#include "BasicGeometry.hpp"
#include "IntervalTree.hpp"
#include "cinolib/meshes/meshes.h"
#include <ctime>
#include <chrono>

namespace fab_translation {
    template <typename T>
    using Vector3 = Eigen::Matrix<T, 3, 1>;

    template <typename T>
    class FabSlicer {
        
    public:
        FabSlicer(mesh::TriMesh<T> tri_mesh, T bottom, T top, T dx,
            T infill_dx)
            : _tri_mesh(tri_mesh), _bottom(bottom), _top(top), _dx(dx), _infill_dx(infill_dx) {

            /* Implement your code here */
            /* 1. Initialize your variables
               2. Build interval tree */
        }

        /* Main entrance for FabSlicer
            return contour and infill_edges, which can be directed sent to corresponding visualization function to visualize in MeshLab
            1. each contour contains a list of point in loop order, each layer can have multiple contours, so the contour is vector<vector<vector<Point>>>.
            2. infill edges is a edge soup for each layer, thus it's vector<vector<Edge>>
        */
        void RunTranslation(std::vector<std::vector<std::vector<Vector3<T>>>>& contour,
            std::vector<std::vector<std::pair<Vector3<T>, Vector3<T>>>>& infill_edges) {
            printf("Start\n");

            std::vector<std::vector<std::pair<indextype, indextype>>> intersection_edges;

            auto t_start = std::chrono::high_resolution_clock::now();
            // Slicing_bruteforce(_tri_mesh, intersection_edges);
            Slicing_accelerated(_tri_mesh, intersection_edges);
            auto t_end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = t_end - t_start;
            printf("Slicing finished: %.6lfs\n", elapsed.count());

            t_start = std::chrono::high_resolution_clock::now();
            CreateContour(_tri_mesh, intersection_edges, contour);
            t_end = std::chrono::high_resolution_clock::now();
            elapsed = t_end - t_start;
            printf("Creating contour finished: %.6lfs\n", elapsed.count());

            t_start = std::chrono::high_resolution_clock::now();
            Infill(contour, infill_edges);
            t_end = std::chrono::high_resolution_clock::now();
            elapsed = t_end - t_start;
            printf("Infilling finished: %.6lfs\n", elapsed.count());
        }

        /* Slicing algorithms
            goal: slice the triangle mesh by a set of parallel planes, 
                  output an intersection edge soup for each layer */
        void Slicing_bruteforce(mesh::TriMesh<T>& tri_mesh, 
            std::vector<std::vector<std::pair<indextype, indextype>>> &intersection_edges) {

            std::vector<Eigen::Vector3i>& elements = tri_mesh.elements();
            std::vector<Vector3<T>>& vertices = tri_mesh.vertices();
            std::vector<Eigen::Vector3i>& edges = tri_mesh.edges();

            intersection_edges.clear();

            for (T h = _bottom; h <= _top; h += _dx) {
                std::vector<std::pair<indextype, indextype>> intersections_one_plane;
                intersections_one_plane.clear();

                geometry::Plane<T> plane(Vector3<T>(0, 0, h), Vector3<T>(0, 0, 1));
                for (int i = 0;i < elements.size();++i) {
                    geometry::Triangle<T> triangle(vertices[elements[i](0)], vertices[elements[i](1)], vertices[elements[i](2)]);
                    triangle.setIndices(elements[i](0), elements[i](1), elements[i](2));
                    std::vector<indextype> intersections = triangle.IntersectPlane(plane);

                    /* Implement your code here */
                    /* What kinds of intersections should be added into intersection edge list? */
                    if ((int)intersections.size() == 2) {
                        intersections_one_plane.push_back(std::make_pair(intersections[0], intersections[1]));
                        intersections_one_plane.push_back(std::make_pair(intersections[1], intersections[0]));
                    }
                }

                intersection_edges.push_back(intersections_one_plane);
            }
        }

        void Slicing_accelerated(mesh::TriMesh<T>& tri_mesh,
            std::vector<std::vector<std::pair<indextype, indextype>>> &intersection_edges) {
            
            std::vector<Eigen::Vector3i>& elements = tri_mesh.elements();
            std::vector<Vector3<T>>& vertices = tri_mesh.vertices();
            std::vector<Eigen::Vector3i>& edges = tri_mesh.edges();

            intersection_edges.clear();

            // for (T h = _bottom; h <= _top; h += _dx) {

            //     std::vector<data_structure::IntervalEntry<T>> candidates;
            //     /* Implement your code here */
            //     /* Retrieve candidate triangle list */

            //     std::vector<std::pair<indextype, indextype>> intersections_one_plane;
            //     intersections_one_plane.clear();

            //     geometry::Plane<T> plane(Vector3<T>(0, 0, h), Vector3<T>(0, 0, 1));
            //     for (int ii = 0;ii < candidates.size();++ii) {
            //         int i = candidates[ii].id;
            //         geometry::Triangle<T> triangle(vertices[elements[i](0)], vertices[elements[i](1)], vertices[elements[i](2)]);
            //         std::vector<indextype> intersections = triangle.IntersectPlane(plane);
                    
            //         /* Implement your code here */
            //         /* What kinds of intersections should be added into intersection edge list? */
            //         if ((int)intersections.size() == 2) {
            //             intersections_one_plane.push_back(std::make_pair(intersections[0], intersections[1]));
            //             intersections_one_plane.push_back(std::make_pair(intersections[1], intersections[0]));
            //         }
            //     }

            //     intersection_edges.push_back(intersections_one_plane);
            // }

            intersection_edges.resize((int)std::floor((_top - _bottom) / _dx));

            // traverse all triangles to calculate intersection
            for (auto &element : elements) {
                T z1 = vertices[element(0)](2);
                T z2 = vertices[element(1)](2);
                T z3 = vertices[element(2)](2);
                T _zmin = std::min(std::min(z1, z2), z3);
                T _zmax = std::max(std::max(z1, z2), z3);
                int idx_bottom = (int)std::ceil((_zmin - _bottom) / _dx);
                int idx_top = (int)std::floor((_zmax - _bottom) / _dx);

                // for each triangle, calculate which set of planes it has intersection with
                geometry::Triangle<T> triangle(vertices[element(0)], vertices[element(1)], vertices[element(2)]);
                triangle.setIndices(element(0), element(1), element(2));
                for (int idx = idx_bottom; idx <= idx_top; idx++) {
                    geometry::Plane<T> plane(Vector3<T>(0, 0, _bottom + _dx * idx), Vector3<T>(0, 0, 1));
                    std::vector<indextype> intersections = triangle.IntersectPlane(plane);
                    if ((int)intersections.size() == 2) {
                        intersection_edges[idx].push_back(std::make_pair(intersections[0], intersections[1]));
                        intersection_edges[idx].push_back(std::make_pair(intersections[1], intersections[0]));
                    }
                }
            }
        }

        /* Find contours
            Goal: Given an intersetion edge soup for each layer, link those edges one by one to form contour loops.
                  Each layer probably has several disjoint contour loops. */
        void CreateContour(mesh::TriMesh<T>& tri_mesh,
            std::vector<std::vector<std::pair<indextype, indextype>>> &intersection_edges,
            std::vector<std::vector<std::vector<Vector3<T>>>>& contours) {
            
            /* Implement your code here */
            /* Input is a edge soup, your task is to generate the loop by linking those edges one by one.
               Thinking about how to find two edge sharing a same end point. set a threshold? or a more clever way? */
            // functions for edge comparison
            std::vector<Vector3<T>>& vertices = tri_mesh.vertices();

            std::vector<std::vector<Vector3<T>>> contours_one_plane;
            std::vector<Vector3<T>> contour;
            T h = _bottom;
            for (auto &edges : intersection_edges) {
                // printf("Current h: %.3lf\n", h);
                contours_one_plane.clear();

                // create forward star representation
                std::map<indextype, std::unordered_set<indextype>> fs;
                fs.clear();
                for (auto &e : edges) {
                    fs[e.first].insert(e.second);
                    fs[e.second].insert(e.first);
                }
                // printf("Points: %d\n", (int)fs.size());

                // find contours
                while (!fs.empty()) {
                    contour.clear();
                    indextype startIdx = fs.begin()->first;
                    indextype curIdx = startIdx, nextIdx;
                    do {
                        // convert index to point coordinates and save them into the contour
                        int idx1 = curIdx >> 32;
                        int idx2 = curIdx & 0xffffffff;
                        // printf("curIdx: %lu, idx1: %d, idx2: %d\n", curIdx, idx1, idx2);
                        if (idx2 == -1)
                            contour.push_back(vertices[idx1]);
                        else {
                            T z1 = vertices[idx1](2), z2 = vertices[idx2](2);
                            T delta = (h - z1) / (z2 - z1);
                            contour.emplace_back(vertices[idx1] + (vertices[idx2] - vertices[idx1]) * delta);
                        }

                        // remove edge pair from forward star data structure
                        auto it = fs.find(curIdx);
                        auto nit = it->second.begin();
                        nextIdx = *nit;
                        it->second.erase(nit);
                        if (it->second.empty())
                            fs.erase(it);
                        auto rit = fs.find(nextIdx);
                        rit->second.erase(rit->second.find(curIdx));
                        if (rit->second.empty())
                            fs.erase(rit);
                        curIdx = nextIdx;
                    }
                    while (curIdx != startIdx);

                    contours_one_plane.push_back(contour);
                    // printf("Found contour with %d points\n", (int)contour.size());
                }

                contours.push_back(contours_one_plane);
                h += _dx;
            }
        }

        /* Generate infill pattern
           Goal: Given the contours at each layer, this function aims to infill the internal part by a pattern which
                 can be procedurally built. (e.g. grid, honey comb, Fermat spiral) 
           The code is for grid pattern, you can rewrite the whole function based on your need. */
        void Infill(std::vector<std::vector<std::vector<Vector3<T>>>>& contours,
            std::vector<std::vector<std::pair<Vector3<T>, Vector3<T>>>>& infill_edges) {
            
            infill_edges.clear();

            for (int i = 0;i < contours.size();++i) {
                std::vector<std::pair<Vector3<T>, Vector3<T>>> infill_edges_one_layer;
                infill_edges_one_layer.clear();

                /* Implement your code here */
                /* 1. find all intersections between contours and your infill pattern 
                2. infill internal space with desired pattern */

                // add infill edges in x direction
                AddInfillEdges(infill_edges_one_layer, contours[i], 0);
                // add infill edges in y direction
                AddInfillEdges(infill_edges_one_layer, contours[i], 1);

                infill_edges.push_back(infill_edges_one_layer);
            }
        }

        // add infill edges along one dimension (x or y)
        // dim - 0: x, 1: y
        void AddInfillEdges(std::vector<std::pair<Vector3<T>, Vector3<T>>> &infill_edges_one_layer,
            std::vector<std::vector<Vector3<T>>> &contours_one_layer, int dim) {

            std::vector<std::vector<Vector3<T>>> intersect_points;

            // determine boundaries
            T _min = 1e10, _max = -1e10;
            for (auto &contour : contours_one_layer)
                for (auto &point : contour) {
                    _min = point(dim) < _min ? point(dim) : _min;
                    _max = point(dim) > _max ? point(dim) : _max;
                }
            int bot = (int)std::ceil(_min / _infill_dx);
            int top = (int)std::floor(_max / _infill_dx);
            intersect_points.resize(top - bot + 1);

            // calculate intersection between contours and infill grid
            for (auto &contour : contours_one_layer) {
                int nPoints = contour.size();
                int idxstart, idxend;
                Vector3<T> p1, p2;
                for (int i = 0, j = 1; i < nPoints; i++, j = (j + 1) % nPoints) {
                    p1 = contour[i], p2 = contour[j];
                    if (p1(dim) > p2(dim))
                        std::swap(p1, p2);
                    if (std::fabs(p2(dim) - p1(dim)) > (T)1e-6) {
                        idxstart = (int)std::ceil(p1(dim) / _infill_dx + (T)1e-6); // leave out p1 to avoid ambiguity
                        idxend = (int)std::floor(p2(dim) / _infill_dx + (T)1e-6);
                        for (int idx = idxstart; idx <= idxend; idx++) {
                            T delta = (_infill_dx * idx - p1(dim)) / (p2(dim) - p1(dim));
                            intersect_points[idx - bot].emplace_back(p1 + (p2 - p1) * delta);
                        }
                    }
                }
            }

            // create infill edges
            for (auto &points : intersect_points) {
                std::sort(points.begin(), points.end(), [=](const auto &p1, const auto &p2) {
                    return p1(dim ^ 1) < p2(dim ^ 1);
                });
                int nPoints = points.size();
                for (int i = 0; i + 1 < nPoints; i += 2)
                    infill_edges_one_layer.push_back(std::make_pair(points[i], points[i + 1]));
            }
        }

        /* Point cloud visualization for each layer's slicing
            file_name: output .ply filename (should be with .ply extension) 
            point_density: the smaller, the denser 
            intersection_edges: edge soup for each layer's slicing */
        void VisualizeSlicing(std::string file_name, 
            T point_density,
            std::vector<std::vector<std::pair<Vector3<T>, Vector3<T>>>> intersection_edges) {
            
            // generate point cloud for ply
            std::vector<Vector3<T>> points;
            points.clear();
            for (int i = 0;i < intersection_edges.size();++i)
                for (int j = 0;j < intersection_edges[i].size();++j) {
                    Vector3<T> s_pos = intersection_edges[i][j].first;
                    Vector3<T> t_pos = intersection_edges[i][j].second;
                    int num_steps = (int)((t_pos - s_pos).norm() / point_density) + 1;
                    for (int step = 0;step <= num_steps;++step) {
                        Vector3<T> pos = s_pos * ((T)step / num_steps) + t_pos * ((T)1.0 - (T)step / num_steps);
                        points.push_back(pos);
                    }
                }

            // output to ply
            FILE* fp = fopen(file_name.c_str(), "w");
            fprintf(fp, "ply\nformat ascii 1.0\n");
            fprintf(fp, "element vertex %d\n", (int)points.size());
            fprintf(fp, "property float32 x\nproperty float32 y\nproperty float32 z\n");
            fprintf(fp, "end_header\n");
            for (int i = 0;i < points.size();++i)
                if (std::is_same<T, float>::value)
                    fprintf(fp, "%.6f %.6f %.6f\n", points[i](0), points[i](1), points[i](2));
                else
                    fprintf(fp, "%.6lf %.6lf %.6lf\n", points[i](0), points[i](1), points[i](2));
            fclose(fp);
        }

        /* Point cloud visualization for each layer's contour
            file_name: output .ply filename (should be with .ply extension) 
            point_density: the smaller, the denser 
            contour: each layer's contour list, each contour is a list a point in loop order */
        void VisualizeContour(std::string file_name,
            T point_density, 
            std::vector<std::vector<std::vector<Vector3<T>>>>& contour) {
            // generate point cloud for ply
            std::vector<Vector3<T>> points;
            points.clear();
            for (int i = 0;i < contour.size();++i)
                for (int j = 0;j < contour[i].size();++j) 
                    for (int k = 0;k < contour[i][j].size();++k) {
                        Vector3<T> s_pos = contour[i][j][k];
                        Vector3<T> t_pos = contour[i][j][(k + 1) % contour[i][j].size()];
                        int num_steps = (int)((t_pos - s_pos).norm() / point_density) + 1;
                        for (int step = 0;step <= num_steps;++step) {
                            Vector3<T> pos = s_pos * ((T)step / num_steps) + t_pos * ((T)1.0 - (T)step / num_steps);
                            points.push_back(pos);
                        }
                    }

            // output to ply
            FILE* fp = fopen(file_name.c_str(), "w");
            fprintf(fp, "ply\nformat ascii 1.0\n");
            fprintf(fp, "element vertex %d\n", (int)points.size());
            fprintf(fp, "property float32 x\nproperty float32 y\nproperty float32 z\n");
            fprintf(fp, "end_header\n");
            for (int i = 0;i < points.size();++i)
                if (std::is_same<T, float>::value)
                    fprintf(fp, "%.6f %.6f %.6f\n", points[i](0), points[i](1), points[i](2));
                else
                    fprintf(fp, "%.6lf %.6lf %.6lf\n", points[i](0), points[i](1), points[i](2));
            fclose(fp);
        }

        /* Point cloud visualization for each layer's slicing
            file_name: output .ply filename (should be with .ply extension) 
            point_density: the smaller, the denser 
            infill_edges: edge soup for each layer's slicing */
        void VisualizeInfill(std::string file_name,
            T point_density,
            std::vector<std::vector<std::pair<Vector3<T>, Vector3<T>>>>& infill_edges) {
            // generate point cloud for ply
            std::vector<Vector3<T>> points;
            points.clear();
            for (int i = 0;i < infill_edges.size();++i)
                for (int j = 0;j < infill_edges[i].size();++j) {
                    int num_steps = (int)((infill_edges[i][j].first - infill_edges[i][j].second).norm() / point_density) + 1;
                    for (int k = 0;k <= num_steps;++k)
                        points.push_back(infill_edges[i][j].first + (infill_edges[i][j].second - infill_edges[i][j].first) * (T)k / (T)num_steps);
                }

            // output to ply
            FILE* fp = fopen(file_name.c_str(), "w");
            fprintf(fp, "ply\nformat ascii 1.0\n");
            fprintf(fp, "element vertex %d\n", (int)points.size());
            fprintf(fp, "property float32 x\nproperty float32 y\nproperty float32 z\n");
            fprintf(fp, "end_header\n");
            for (int i = 0;i < points.size();++i)
                if (std::is_same<T, float>::value)
                    fprintf(fp, "%.6f %.6f %.6f\n", points[i](0), points[i](1), points[i](2));
                else
                    fprintf(fp, "%.6lf %.6lf %.6lf\n", points[i](0), points[i](1), points[i](2));
            fclose(fp);
        }

    private:
        mesh::TriMesh<T> _tri_mesh;

        /* Variables for slicing */
        T _bottom, _top, _dx;

        /* Variables for infill algorithm */
        T _infill_dx;                                   // infill pattern will be equal-length grid
        T _infill_x_lower_bound, _infill_x_upper_bound;
        T _infill_y_lower_bound, _infill_y_upper_bound;

        /* accelerated data structure */
        data_structure::IntervalTree<T> _interval_tree;

    };
}