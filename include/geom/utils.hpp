//
// Created by xetql on 3/17/21.
//

#pragma once
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Point_2.h>
#include <CGAL/Vector_2.h>
#include <CGAL/Ray_2.h>
#include <CGAL/iterator.h>
#include <CGAL/ch_jarvis.h>
// Kernels
#include <CGAL/Cartesian.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/ch_graham_andrew.h>
#include <CGAL/Polygon_2_algorithms.h>

#include <vector>
#include <array>
#include <cassert>
#include "algorithm.hpp"
#include "numeric/utils.hpp"

//using K      = CGAL::Simple_cartesian<double>;
using K        = CGAL::Exact_predicates_inexact_constructions_kernel;
using ExactK   = CGAL::Exact_predicates_exact_constructions_kernel;
using Polygon2 = CGAL::Polygon_2<K>;
using Line2    = CGAL::Line_2<K>;
using Segment2 = CGAL::Segment_2<K>;
using Point2   = CGAL::Point_2<K>;
using EPoint2  = CGAL::Point_2<ExactK>;
using Vector2  = CGAL::Vector_2<K>;
using Ray2     = CGAL::Ray_2<K>;
template<class Real>
Real side(Real tx, Real ty, Real dx, Real dy) {
    return dx * ty - dy * tx;
}
template<class Real>
void add_to_bisection(std::vector<Point2> &b1, std::vector<Point2> &b2, Real vx, Real vy, const Point2 &pmed, const Point2 &p) {
    auto median_to_p = p-pmed;
    auto s = sign(CGAL::to_double(side(vx, vy, (Real) median_to_p.x(), (Real) median_to_p.y())));
    if (s <= 0) {
        b1.push_back(p);
    } else {
        b2.push_back(p);
    }
}
inline auto get_angle(const Vector2 &v, const Vector2 &origin) {
    auto dot = v * origin;
    auto det = v.x() * origin.y() - v.y() * origin.x();
    return std::atan2(CGAL::to_double(det), CGAL::to_double(dot));
}

template<class Real>
inline auto get_angle(Real x1, Real y1, Real x2, Real y2) {
    auto dot = x1*x2 + y1*y2;
    auto det = x1 * y2 - y1 * x2;
    return std::atan2(det, dot);
}

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

template<class Real>
constexpr Real to_degree(Real radian){
    // 57.29577951308232 = 180 / pi
    return 57.29577951308232 * radian;
}

constexpr auto to_radian(double degree){
    // 0.017453292519943295 = PI / 180
    return 0.017453292519943295 * degree;
}

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
        return almost_equal(CGAL::to_double(p1.x()), CGAL::to_double(p2.x()), 2) &&
        almost_equal(CGAL::to_double(p1.y()), CGAL::to_double(p2.y()), 2);
    }
};
template<class Real>
std::pair<Polygon2, Polygon2> bisect_polygon(const Polygon2 &poly, Real vx, Real vy, const Point2 &median) {

	CGAL::Vector_2<ExactK> evec(vx, vy);
    
	CGAL::Point_2<ExactK > emedian(CGAL::to_double(median.x()), CGAL::to_double(median.y()));

    CGAL::Line_2<ExactK> med(emedian, evec);

    std::vector<Point2> intersections {};
    /* iterate over polygon edges to find intersections with separating line with exact predicates...*/
    for (auto eit = poly.edges_begin(); eit != poly.edges_end(); ++eit) {
        auto nonExactSegment = *eit;
        CGAL::Segment_2<ExactK> s(EPoint2(nonExactSegment.source().x(), nonExactSegment.source().y()),
                                  EPoint2(nonExactSegment.target().x(), nonExactSegment.target().y()));

        // Compute intersection with *positive* ray
    
		auto inter = CGAL::intersection(s, med);

        // if we find an intersection, then add the intersection point to the list of intersections
        if (inter.has_value()) {
            auto ep = boost::get<EPoint2>(inter.value());
            Point2 p(CGAL::to_double(ep.x()), CGAL::to_double(ep.y()));
            intersections.push_back(p);
        }
    }

    // remove duplicates
    auto last = distinct(intersections.begin(), intersections.end(), P2Comp{});
    intersections.erase(last, intersections.end());

    if(intersections.size() != 2) {
        throw std::logic_error("A line intersects 2 polygon edges.");
    }

    std::vector<Point2> b1(intersections.begin(), intersections.end()),
    b2(intersections.begin(), intersections.end());

    // add to left or right
    for (auto eit = poly.edges_begin(); eit != poly.edges_end(); ++eit) {
        Segment2 s = *eit;
        add_to_bisection(b1, b2, vx, vy, median, s.source());
        add_to_bisection(b1, b2, vx, vy, median, s.target());
    }

    std::vector<Point2> chull1{}, chull2{};

    CGAL::ch_jarvis(b1.begin(), b1.end(), std::back_inserter(chull1));
    CGAL::ch_jarvis(b2.begin(), b2.end(), std::back_inserter(chull2));

    Polygon2 poly1(chull1.begin(), chull1.end());

    Polygon2 poly2(chull2.begin(), chull2.end());

    return {poly1, poly2};
}
