#pragma once
#include <Eigen/Dense>
#include <vector>
#include <iostream>

// index type
typedef unsigned long long indextype;

namespace geometry {
    template <typename T>
    using Vector3 = Eigen::Matrix<T, 3, 1>;

    template <typename T>
    using Matrix3 = Eigen::Matrix<T, 3, 3>;

    // the plane is represented by (x - _p) /dot _normal = 0
    template <typename T>
    class Plane {
    public:
        Plane(Vector3<T> p, Vector3<T> normal) {
            _p = p;
            _normal = normal;
            _normal.normalize();
        }

        Vector3<T>& p() { return _p; }
        Vector3<T>& normal() { return _normal; }
        
        // return if the point is on plane
        // also fill parameter dist as the signed distance from point to plane
        bool onPlane(Vector3<T> point, T& dist) {
            dist = (point - _p).dot(_normal);
            if (std::fabs(dist) < 1e-6) {
                return true;
            } else {
                return false;
            }
        }

    private:
        Vector3<T> _p;
        Vector3<T> _normal;
    };

    template <typename T>
    class Triangle {
    public:
        Triangle(Vector3<T> v0, Vector3<T> v1, Vector3<T> v2) {
            // Pre-order vertices according to Y coordinate
            if (v0(0) > v1(0)) std::swap(v0, v1);
            if (v0(0) > v2(0)) std::swap(v0, v2);
            if (v1(0) > v2(0)) std::swap(v1, v2);

            _vertices[0] = v0;
            _vertices[1] = v1;
            _vertices[2] = v2;

            // Add normal calculation
            updateDeltas();
            updateNormal();
            updateBoundaries();
        }

        Vector3<T>* vertices() { return _vertices; }
        Vector3<T>& vertices(int idx) { return _vertices[idx]; }

        void setIndices(int i1, int i2, int i3) {
            _idx[0] = i1;
            _idx[1] = i2;
            _idx[2] = i3;
        }
        int *indices() { return _idx; }
        int &indices(int idx) { return _idx[idx]; }

        inline Vector3<T> normal() { return _normal; }
        inline void updateNormal() { _normal = std::move(-_delta[0].cross(_delta[2])); }

        inline Vector3<T> delta(int idx) { return _delta[idx]; }
        inline void updateDeltas() {
            _delta[0] = std::move(_vertices[1] - _vertices[0]);
            _delta[1] = std::move(_vertices[2] - _vertices[1]);
            _delta[2] = std::move(_vertices[0] - _vertices[2]);
        }

        inline Vector3<T> bmin() { return _bmin; }
        inline Vector3<T> bmax() { return _bmax; }
        inline void updateBoundaries() {
            _bmin = std::move(_vertices[0].cwiseMin(_vertices[1]).cwiseMin(_vertices[2]));
            _bmax = std::move(_vertices[0].cwiseMax(_vertices[1]).cwiseMax(_vertices[2]));
        }

        /* Implement triangle plane intersection.
            Input is a plane, output is a list of intersection points
            Hints:
                1. Take some considerations: handle duplicate intersections (when plane intersect with 
                    triangle at corner but you find them by edge-plane intersection).
                2. You can modify the function input/return type to any format you want.
        */
        std::vector<indextype> IntersectPlane(Plane<T> p) {
            /* Implement your code here */
            T dist[3];
            for (int i = 0; i < 3; i++)
                dist[i] = ftz((_vertices[i] - p.p()).dot(p.normal()));

            // compute intersection points
            std::vector<indextype> points;
            points.clear();
            for (int i = 0, j = 1; i < 3; i++, j = (j + 1) % 3) {
                if (dist[i] == 0.0)
                    points.push_back((indextype)_idx[i] << 32 ^ 0xffffffff);
                if (dist[j] == 0.0)
                    points.push_back((indextype)_idx[j] << 32 ^ 0xffffffff);
                if (dist[i] * dist[j] < 0.0)
                    if (_idx[i] < _idx[j])
                        points.push_back((indextype)_idx[i] << 32 ^ _idx[j]);
                    else
                        points.push_back((indextype)_idx[j] << 32 ^ _idx[i]);
            }

            // remove duplicate points
            std::vector<indextype> result;
            result.clear(); 
            int nPoints = points.size();
            for (int i = 0, j = 1; i < nPoints; i++, j = (j + 1) % nPoints)
                if (points[i] != points[j])
                    result.push_back(points[i]);
            
            return std::move(result);
        }

        // Assignment 2: implement ray-triangle intersection.
        // The ray is defined as r(t) = origin + t * dir.
        // You should return a scalar t such that r(t) is the intersection point. Which value
        // to return for the case of no intersections is up to you. You can also change the
        // signature of this function as you see fit.
        const T IntersectRay(const Vector3<T>& origin, const Vector3<T>& dir) const {
            /* Assignment 2. */
            /* Implement your code here */
            T t = -(origin - _vertices[0]).dot(_normal) / dir.dot(_normal);
            if (t < 0) return -1e10;

            // Check whether intersection point is located within the triangle or not
            Vector3<T> target = std::move(origin + t * dir);
            T d1 = (target - _vertices[0]).cross(_delta[0]).dot(_normal);
            T d2 = (target - _vertices[1]).cross(_delta[1]).dot(_normal);
            T d3 = (target - _vertices[2]).cross(_delta[2]).dot(_normal);
            return (d1 * d2 >= 0.0 && d1 * d3 >= 0.0 && d2 * d3 >= 0.0) ? t : -1e10;
        }

        // Rotate the triangle by quaternion q
        void rotate(const Eigen::Quaternion<T> &q) {
            Matrix3<T> mat = std::move(q.toRotationMatrix());
            for (int i = 0; i < 3; ++i)
                _vertices[i] = mat * _vertices[i];

            if (_vertices[0](0) > _vertices[1](0)) std::swap(_vertices[0], _vertices[1]);
            if (_vertices[0](0) > _vertices[2](0)) std::swap(_vertices[0], _vertices[2]);
            if (_vertices[1](0) > _vertices[2](0)) std::swap(_vertices[1], _vertices[2]);

            updateDeltas();
            updateNormal();
            updateBoundaries();
        }
        
    private:
        inline T ftz(T x) { return std::fabs(x) < 1e-6 ? 0 : x; }

        Vector3<T> _vertices[3];
        Vector3<T> _delta[3];
        int _idx[3];

        // Norm of triangle
        Vector3<T> _normal;

        // Bounding box
        Vector3<T> _bmin, _bmax;
    };
}
