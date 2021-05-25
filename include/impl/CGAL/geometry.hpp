//
// Created by xetql on 5/25/21.
//

#ifndef YALBB_GEOMETRY_HPP
#define YALBB_GEOMETRY_HPP
#include "geometry_utils.hpp"
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Point_2.h>
#include <CGAL/Vector_2.h>
#include <CGAL/Ray_2.h>
#include <CGAL/iterator.h>

// Kernels
#include <CGAL/Cartesian.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/ch_graham_andrew.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/ch_jarvis.h>

using K        = CGAL::Exact_predicates_inexact_constructions_kernel;
using ExactK   = CGAL::Exact_predicates_exact_constructions_kernel;
using Polygon2 = CGAL::Polygon_2<K>;
using Segment2 = CGAL::Segment_2<K>;
using Point2   = CGAL::Point_2<K>;
using EPoint2  = CGAL::Point_2<ExactK>;
using Vector2  = CGAL::Vector_2<K>;

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

inline auto side(const Vector2 &T, const Vector2 &D) -> std::remove_const_t<std::remove_reference_t<decltype(T.y())>> {
    return D.x() * T.y() - D.y() * T.x();
}

inline void add_to_bisection(std::vector<Point2> &b1, std::vector<Point2> &b2, const Vector2 &v, const Point2 &pmed, const Point2 &p) {
    auto s = sign(CGAL::to_double(side(v, p - pmed)));
    if (s <= 0) {
        b1.push_back(p);
    } else {
        b2.push_back(p);
    }
}

struct P2Comp {
    int operator()(Point2& p1, Point2& p2) {
        return almost_equal(CGAL::to_double(p1.x()), CGAL::to_double(p2.x()), 2) &&
               almost_equal(CGAL::to_double(p1.y()), CGAL::to_double(p2.y()), 2);
    }
};

template<class Real>
Polygon2 init_domain(Real minx, Real miny, Real maxx, Real maxy) {
    Point2 p1(minx, miny), p2(maxx, miny), p3(maxx, maxy), p4(minx, maxy);
    Polygon2 d;
    d.push_back(p1);
    d.push_back(p2);
    d.push_back(p3);
    d.push_back(p4);
    return d;
}


inline std::pair<Polygon2, Polygon2> bisect_polygon(const Polygon2 &poly, const Vector2& vec, const Point2 &median) {

    CGAL::Vector_2<ExactK> evec(CGAL::to_double(vec.x()), CGAL::to_double(vec.y()));

    CGAL::Point_2<ExactK > emedian(CGAL::to_double(median.x()), CGAL::to_double(median.y()));

    CGAL::Line_2<ExactK> med(emedian, evec);

    std::vector<Point2> intersections {};
    // iterate over polygon edges to find intersections with separating line with exact predicates...
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
    // iterate over polygon edges to find intersections with separating line with exact predicates...
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

#endif //YALBB_GEOMETRY_HPP
