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
#include <serial/format.hpp>
#include <ostream>
#include <parallel/algorithm.hpp>
#include "algorithm.hpp"
#include "numeric/utils.hpp"


using K        = CGAL::Exact_predicates_inexact_constructions_kernel;
using ExactK   = CGAL::Exact_predicates_exact_constructions_kernel;

using Polygon2 = CGAL::Polygon_2<K>;

using Segment2 = CGAL::Segment_2<K>;
using Point2   = CGAL::Point_2<K>;
using EPoint2  = CGAL::Point_2<ExactK>;
using Vector2  = CGAL::Vector_2<K>;

static const double acceptable_error = std::numeric_limits<double>::epsilon() * 16;

struct Point {
    long double x, y;

    Point operator-(const Point& p) const {
        return Point {x-p.x, y-p.y};
    }

    Point operator+(const Point& p) const {
        return Point {x+p.x, y+p.y};
    }

    friend std::ostream &operator<<(std::ostream &os, const Point &point) {
        os << "x: " << point.x << " y: " << point.y;
        return os;
    }

    bool operator==(const Point &rhs) const {
        return almost_equal(x, rhs.x, 4) && almost_equal(y, rhs.y, 4);
    }

    bool operator!=(const Point &rhs) const {
        return !(rhs == *this);
    }
};
using Vector = Point;

struct Segment {
    Point a, b;

    Segment operator-(const Point &p) const {
        return Segment { Point {a.x - p.x, a.y - p.y}, Point {b.x - p.x, b.y - p.y} };
    }

    [[nodiscard]] bool contains (const Point& p) const {

        const auto crossproduct = (p.y - a.y) * (b.x - a.x) - (p.x - a.x) * (b.y - a.y);

        if(std::fabs(crossproduct) > acceptable_error) {
            return false;
        }

        const auto dotproduct = (p.x - a.x) * (b.x - a.x) + (p.y - a.y) * (b.y - a.y);

        if(dotproduct < 0.0) {
            return false;
        }

        const auto squaredlengthab = (b.x - a.x) * (b.x - a.x) + (b.y - a.y) * (b.y - a.y);

        if(dotproduct > squaredlengthab) {
            return false;
        }

        return true;
    }

    friend std::ostream &operator<<(std::ostream &os, const Segment &segment) {
        os << "a: " << segment.a << " b: " << segment.b;
        return os;
    }
};

struct Line {
    Point p;
    Vector v;
    explicit Line(const Segment& s) : p(s.b), v(s.a - s.b) {}
    Line(Point&& p, Vector&& v) : p(p), v(v) {}
    Line(const Point& p, const Vector& v) : p(p), v(v) {}

    [[nodiscard]] std::pair<Point, Point> get_points() const {
        return {p, p+v};
    }
};

inline std::optional<Point> intersection(const Line& l1, const Line& l2) {
    const auto&[p1, p2] = l1.get_points();
    const auto&[p3, p4] = l2.get_points();

    const auto D = (p1.x-p2.x)*(p3.y-p4.y) - (p1.y-p2.y)*(p3.x-p4.x);

    if(std::fabs(D) < acceptable_error) {
        return std::nullopt;
    }

    const auto x = ((p1.x*p2.y-p1.y*p2.x)*(p3.x-p4.x) - (p1.x-p2.x)*(p3.x*p4.y - p3.y*p4.x)) / D;
    const auto y = ((p1.x*p2.y-p1.y*p2.x)*(p3.y-p4.y) - (p1.y-p2.y)*(p3.x*p4.y - p3.y*p4.x)) / D;

    return Point{x, y};
}

inline std::optional<Point> intersection(const Line& l, const Segment& s) {
    const Line l_s(s);
    const auto opt_i = intersection(l, l_s);
    if(opt_i.has_value()){
        const auto& i = opt_i.value();
        if(s.contains(i))
            return i;
        else {
            return std::nullopt;
        }
    } else {
        return std::nullopt;
    }
}

template<class Real>
Real side(Real tx, Real ty, Real dx, Real dy) {
    return dx * ty - dy * tx;
}

inline void add_to_bisection(std::vector<Point> &b1, std::vector<Point> &b2, const Vector& v, const Point &pmed, const Point &p) {
    const auto median_to_p = p - pmed;
    const auto vo = v-pmed;
    const auto s = sign(side(v.x, v.y, median_to_p.x, median_to_p.y));

    if (s <= 0) {
        b1.push_back(p);
    } else {
        b2.push_back(p);
    }
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
    int operator()(Point2& p1, Point2& p2) {
        return almost_equal(CGAL::to_double(p1.x()), CGAL::to_double(p2.x()), 2) &&
        almost_equal(CGAL::to_double(p1.y()), CGAL::to_double(p2.y()), 2);
    }
};

inline std::pair<Polygon2, Polygon2> bisect_polygon(const Polygon2 &poly, const Vector2& vec, const Point2 &median) {

    CGAL::Vector_2<ExactK> evec(CGAL::to_double(vec.x()), CGAL::to_double(vec.y()));

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
        add_to_bisection(b1, b2, vec, median, s.source());
        add_to_bisection(b1, b2, vec, median, s.target());
    }

    std::vector<Point2> chull1{}, chull2{};

    CGAL::ch_jarvis(b1.begin(), b1.end(), std::back_inserter(chull1));
    CGAL::ch_jarvis(b2.begin(), b2.end(), std::back_inserter(chull2));

    Polygon2 poly1(chull1.begin(), chull1.end());

    Polygon2 poly2(chull2.begin(), chull2.end());

    return {poly1, poly2};
}

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

template<class Real>
std::pair<Polygon2, Polygon2> bisect_polygon(const Polygon2 &poly, Real vx, Real vy, Real px, Real py) {

    Vector vec {vx, vy};
    Point  pmed{px, py};
    Line   med(pmed, vec);

    std::vector<Point> intersections {};

    /* iterate over polygon edges to find intersections with separating line...*/
    for (auto eit = poly.edges_begin(); eit != poly.edges_end(); ++eit) {
        auto cgal_segment = *eit;
        const Segment s {Point{cgal_segment.source().x(), cgal_segment.source().y()},
                   Point{cgal_segment.target().x(), cgal_segment.target().y()}
        };

        // Compute intersection
        const auto opt = intersection(med, s);

        // if we find an intersection, then add the intersection point to the list of intersections
        if (opt.has_value()) {
            const auto& inter = opt.value();
            intersections.push_back(inter);
        }
    }

    // remove duplicates
    /*if(intersections.size() > 2){
        auto last = distinct(intersections.begin(), intersections.end(), std::equal_to{});
        intersections.erase(last, intersections.end());
    }*/

    if (intersections.size() != 2) {
        std::stringstream ss {};
        std::for_each(intersections.begin(), intersections.end(), [&ss](auto p) {ss << p;});
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if(!rank) throw std::logic_error(ser::fmt("A line intersects 2 polygon edges. %s", ss.str()));
    }

    std::vector<Point> b1(intersections.begin(), intersections.end()), b2(intersections.begin(), intersections.end());
    b1.reserve(b1.size() + poly.size());
    b2.reserve(b2.size() + poly.size());

    // add to left or right
    for (auto eit = poly.edges_begin(); eit != poly.edges_end(); ++eit) {
        const Segment2 s = *eit;
        add_to_bisection(b1, b2, vec, pmed, Point{s.source().x(), s.source().y()});
        add_to_bisection(b1, b2, vec, pmed, Point{s.target().x(), s.target().y()});
    }

    std::vector<Point2> b1_cgal{}, b2_cgal{};
    b1_cgal.reserve(b1.size());
    b2_cgal.reserve(b2.size());

    std::transform(b1.begin(), b1.end(), std::back_inserter(b1_cgal), [](const auto& p){
        return Point2(p.x, p.y);
        //return Point2(p.x.template convert_to<double>(), p.y.template convert_to<double>());
    });
    std::transform(b2.begin(), b2.end(), std::back_inserter(b2_cgal), [](const auto& p){
        return Point2(p.x, p.y);
        //return Point2(p.x.template convert_to<double>(), p.y.template convert_to<double>());
    });

    std::vector<Point2> chull1{}, chull2{};

    CGAL::ch_jarvis(b1_cgal.begin(), b1_cgal.end(), std::back_inserter(chull1));
    CGAL::ch_jarvis(b2_cgal.begin(), b2_cgal.end(), std::back_inserter(chull2));

    Polygon2 poly1(chull1.begin(), chull1.end());
    Polygon2 poly2(chull2.begin(), chull2.end());

    return {poly1, poly2};
}
