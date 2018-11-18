#pragma once
#include <Eigen/Dense>
#include <vector>
#include <set>
#include <iostream>
#include <functional>

namespace geometry {
    template <typename T>
    using Vector3 = Eigen::Matrix<T, 3, 1>;
    template <typename T>
    using Vector2 = Eigen::Matrix<T, 2, 1>;
    template <typename T>
    using VectorX = Eigen::Matrix<T, Eigen::Dynamic, 1>;

    template <typename T>
    class VectorCompare {
    public:
        bool operator () (const Vector3<T> &v1, const Vector3<T> &v2) { return v1(1) < v2(1); }
    };

    // Q1: implement Graham Scan for 2D convex hull
    // The input is a vector of 2D points
    // The output should be the 2D points on the convex hull
    // Remove duplicates and sort your return vector
    // for comparing with the reference result
    template <typename T>
    std::vector<Vector2<T>> ConvexHull2D(const std::vector<Vector2<T>> &points) {
        std::vector<Vector2<T>> vec(points);
        auto it = std::min_element(vec.begin(), vec.end(),
            [](const Vector2<T> &v1, const Vector2<T> &v2) { return v1(1) < v2(1); });

        // sorting by intersection angles
        Vector2<T> V0 = *it;
        vec.erase(it);
        std::sort(vec.begin(), vec.end(), [&](const Vector2<T> &v1, const Vector2<T> &v2) {
            double cos1 = (v1(0) - V0(0)) / (v1 - V0).norm();
            double cos2 = (v2(0) - V0(0)) / (v2 - V0).norm();
            return cos1 < cos2 || cos1 == cos2 && (v1 - V0).norm() < (v2 - V0).norm();
        });
        vec.push_back(V0);

        // computing convex hull
        std::vector<Vector2<T>> C;
        C.push_back(V0);
        C.push_back(vec[0]);
        auto cross = [](const Vector2<T> &v1, const Vector2<T> &v2) {
            return v1(0) * v2(1) - v1(1) * v2(0);
        };
        for (auto &v : vec) {
            while ((int)C.size() > 1 && cross(v - C.back(), C.back() - C[C.size() - 2]) <= 0)
                C.pop_back();
            C.push_back(v);
        }

        // remove duplicates and sort result
        C.pop_back();
        std::sort(C.begin(), C.end(), [](const Vector2<T> &v1, const Vector2<T> &v2) {
            return v1(0) < v2(0) || v1(0) == v2(0) && v1(1) < v2(1);
        });
        return std::move(C);
    }

    // Q2: implement brute-force method for Nd Pareto front
    // The input is a vector of Nd points
    // The output should be the Nd points on the Pareto front
    // sort your return vector for comparing with 
    // the reference result
    template <typename T>
    std::vector<VectorX<T>> ParetoFrontNdNaive(const std::vector<VectorX<T>> &points) {
        const int N = points.size();
        auto vectorLE = [](const VectorX<T> &v1, const VectorX<T> &v2) {
            bool ret = true;
            for (int i = 0; i < v1.rows(); ++i)
                ret &= v1(i) <= v2(i);
            return ret;
        };

        std::vector<VectorX<T>> P;
        for (int i = 0; i < N; ++i) {
            bool flag = true;
            for (int j = 0; j < N; ++j)
                flag &= i == j || !vectorLE(points[j], points[i]);
            if (flag)
                P.push_back(points[i]);
        }

        std::function<bool(const VectorX<T> &, const VectorX<T> &)> compare;
        compare = [&](const VectorX<T> &v1, const VectorX<T> &v2) -> bool {
            if (v1(0) < v2(0))
                return true;
            else if (v1(0) > v2(0))
                return false;
            else
                return v1.rows() > 1 ? compare(v1.tail(v1.rows() - 1), v2.tail(v2.rows() - 1)) : false;
        };
        std::sort(P.begin(), P.end(), compare);

        return std::move(P);
    }

    // Q3: implement nlog(n) method for 2d Pareto front
    // The input is a vector of 2d points
    // The output should be the 2d points on the Pareto front
    // sort your return vector for comparing with 
    // the reference result
    template <typename T>
    std::vector<Vector2<T>> ParetoFront2D(const std::vector<Vector2<T>> &points) {
        std::vector<Vector2<T>> vec(points);
        std::sort(vec.begin(), vec.end(), [](const Vector2<T> &v1, const Vector2<T> &v2) {
            return v1(0) < v2(0) || v1(0) == v2(0) && v1(1) < v2(1);
        });

        std::vector<Vector2<T>> P;
        P.push_back(vec[0]);
        for (auto &v : vec)
            if (v(1) < P.back()(1))
                P.push_back(v);

        return std::move(P);
    }

    // bonus question 1: implement 3D convex hull
    // The input is a vector of 3D points
    // The output should be the 3D points on the convex hull
    // sort your return vector for comparing with 
    // the reference result
    template <typename T>
    std::vector<Vector3<T>> ConvexHull3D(const std::vector<Vector3<T>> &points) {
        return points;
    }

    // bonus question 2: implement nlog(n) method for 3d Pareto front
    // The input is a vector of 3d points
    // The output should be the 3d points on the Pareto front
    // sort your return vector for comparing with 
    // the reference result
    template <typename T>
    std::vector<Vector3<T>> ParetoFront3D(const std::vector<Vector3<T>> &points) {
        std::vector<Vector3<T>> vec(points);
        std::sort(vec.begin(), vec.end(), [](const Vector3<T> &a, const Vector3<T> &b) {
            return a(0) < b(0) || a(0) == b(0) && (a(1) < b(1) || a(1) == b(1) && a(2) < b(2));
        });

        std::set<Vector3<T>, VectorCompare<T>> front;
        std::vector<Vector3<T>> P;
        for (auto &vi : vec) {
            auto it = std::lower_bound(front.begin(), front.end(), vi, VectorCompare<T>());
            if ((it == front.begin() || vi(2) < (*std::prev(it, 1))(2)) &&
                (it == front.end() || vi(1) < (*it)(1) || vi(1) == (*it)(1) && vi(2) < (*it)(2))) {
                P.push_back(vi);
                while (it != front.end() && vi(1) <= (*it)(1) && vi(2) <= (*it)(2))
                    it = front.erase(it);
                front.insert(it == front.begin() ? it : std::prev(it, 1), vi);
            }
        }

        return std::move(P);
    }
}
