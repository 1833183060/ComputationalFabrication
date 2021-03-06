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
        // also fill parameter dist field as the signed distance from point to plane
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
            _vertices[0] = v0;
            _vertices[1] = v1;
            _vertices[2] = v2;

            // Add normal calculation
            updateNormal();
            // calcInvMatrix();
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

        Vector3<T> normal() { return _normal; }
        void updateNormal() { _normal = (_vertices[1] - _vertices[0]).cross(_vertices[2] - _vertices[0]); }
        void updateInvMatrix() { calcInvMatrix(); }

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
            Vector3<T> target = origin + t * dir;
            T d1 = (target - _vertices[0]).cross(_vertices[1] - _vertices[0]).dot(_normal);
            T d2 = (target - _vertices[1]).cross(_vertices[2] - _vertices[1]).dot(_normal);
            T d3 = (target - _vertices[2]).cross(_vertices[0] - _vertices[2]).dot(_normal);
            return (d1 * d2 >= 0.0 && d1 * d3 >= 0.0 && d2 * d3 >= 0.0) ? t : -1e10;
        }
        
    private:
        inline T ftz(T x) { return std::fabs(x) < 1e-6 ? 0 : x; }

        inline void calcInvMatrix() {
            Matrix3<T> _coords;
            for (int i = 0; i < 3; i++)
                _coords.col(i) = _vertices[i];
            _inv = _coords.inverse();
        }

        Vector3<T> _vertices[3];
        int _idx[3];

        // Norm of triangle
        Vector3<T> _normal;

        // Inverted matrix
        Matrix3<T> _inv;
    };
}
