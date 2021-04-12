//
// Created by xetql on 3/17/21.
//

#pragma once
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Point_2.h>
#include <CGAL/Vector_2.h>
#include <CGAL/Ray_2.h>
#include <CGAL/iterator.h>
#include <CGAL/ch_graham_andrew.h>
#include <CGAL/Polygon_2_algorithms.h>

#include <vector>
#include <array>

#include "numeric/utils.hpp"

using K        = CGAL::Exact_predicates_inexact_constructions_kernel;
using Polygon2 = CGAL::Polygon_2<K>;
using Segment2 = CGAL::Segment_2<K>;
using Point2   = CGAL::Point_2<K>;
using Vector2  = CGAL::Vector_2<K>;
using Ray2     = CGAL::Ray_2<K>;

inline auto get_angle(const Vector2 &v, const Vector2 &origin) {
    auto dot = v * origin;
    auto det = v.x() * origin.y() - v.y() * origin.x();
    return std::atan2(det, dot);
}

std::pair<Polygon2, Polygon2> bisect_polygon(const Polygon2& poly, const Vector2& vec, const Point2& median);

void add_to_bisection(std::vector<Point2> &b1, std::vector<Point2> &b2, const Vector2 &v, const Point2 &pmed, const Point2 &p);

template<class Real>
void rotate(const std::array<std::array<Real, 2>, 2> &matrix, std::vector<Real> &x, std::vector<Real> &y){
    const auto size = x.size();
    assert(x.size() == y.size());
    Real rxi, ryi;
#pragma GCC ivdep
    for (unsigned i = 0; i < size; ++i) {
        rxi = matrix[0][0] * x[i] + matrix[0][1] * y[i];
        ryi = matrix[1][0] * x[i] + matrix[1][1] * y[i];
        x[i] = rxi;
        y[i] = ryi;
    }
}
template<class Real>
std::pair<Real, Real> rotate(const std::array<std::array<Real, 2>, 2> &matrix, Real x, Real y) {
    return {matrix[0][0] * x + matrix[0][1] * y, matrix[1][0] * x + matrix[1][1] * y};
}

std::pair<std::vector<double>, std::vector<double>> rotate_copy(
        const std::array<std::array<double, 2>, 2> &matrix,
        const std::vector<double> &x,
        const std::vector<double> &y);

auto side(const Vector2 &T, const Vector2 &D) -> std::remove_const_t<std::remove_reference_t<decltype(T.y())>>;

template<class Real> std::array<std::array<Real, 2>, 2> get_rotation_matrix(Real theta) {
    std::array<std::array<Real, 2>, 2> rotate{};
    rotate[0][0] =  std::cos(theta);
    rotate[0][1] = -std::sin(theta);
    rotate[1][0] =  std::sin(theta);
    rotate[1][1] =  std::cos(theta);
    return rotate;
}

template<class ForwardIt>
void rotate(const std::array<std::array<typename ForwardIt::value_type, 2>, 2> &matrix, ForwardIt x_begin, ForwardIt x_end, ForwardIt y_begin) {
    const auto size = std::distance(x_begin, x_end);
    typename ForwardIt::value_type rxi, ryi;
#pragma GCC ivdep
    for (auto i = 0; i < size; ++i) {
        rxi = matrix[0][0] * (*(x_begin + i)) + matrix[0][1] * (*(y_begin + i));
        ryi = matrix[1][0] * (*(x_begin + i)) + matrix[1][1] * (*(y_begin + i));
        (*(x_begin + i)) = rxi;
        (*(y_begin + i)) = ryi;
    }
}
template<class Real, class ForwardIt, class GetPositionFunc>
void rotate(const std::array<std::array<Real, 2>, 2> &rotation_mat, ForwardIt el_begin, ForwardIt el_end, GetPositionFunc getPosition) {
    const auto size = std::distance(el_begin, el_end);
    Real rxi, ryi;
#pragma GCC ivdep
    for (unsigned i = 0; i < size; ++i) {
        auto position = getPosition(&(*(el_begin+i)));
        rxi = rotation_mat[0][0] * (position->at(0)) + rotation_mat[0][1] * (position->at(1));
        ryi = rotation_mat[1][0] * (position->at(0)) + rotation_mat[1][1] * (position->at(1));
        (position->at(0)) = rxi;
        (position->at(1)) = ryi;
    }
}

struct P2Comp {
    int operator()(Point2 p1, Point2 p2) {
        return almost_equal(p1.x(), p2.x(), 2) && almost_equal(p1.y(), p2.y(), 2);
    }
};