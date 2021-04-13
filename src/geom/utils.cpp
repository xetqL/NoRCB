//
// Created by xetql on 3/17/21.
//
#include "geom/utils.hpp"
#include "algorithm.hpp"

auto side(const Vector2 &T, const Vector2 &D) -> std::remove_const_t<std::remove_reference_t<decltype(T.y())>> {
    return D.x() * T.y() - D.y() * T.x();
}

std::pair<std::vector<double>, std::vector<double>> rotate_copy(const std::array<std::array<double, 2>, 2> &matrix, const std::vector<double> &x,
                                                                const std::vector<double> &y) {
    const auto size = x.size();
    std::vector<double> rx(size, 0), ry(size, 0);
#pragma GCC ivdep
    for (unsigned i = 0; i < size; ++i) {
        rx[i] = matrix[0][0] * x[i] + matrix[0][1] * y[i];
        ry[i] = matrix[1][0] * x[i] + matrix[1][1] * y[i];
    }
    return {rx, ry};
}

void add_to_bisection(std::vector<Point2> &b1, std::vector<Point2> &b2, const Vector2 &v, const Point2 &pmed, const Point2 &p) {
    auto s = sign(side(v, p - pmed));
    if (s <= 0) {
        b1.push_back(p);
    } else {
        b2.push_back(p);
    }
}

std::pair<Polygon2, Polygon2> bisect_polygon(const Polygon2 &poly, const Vector2 &vec, const Point2 &median) {
    Ray2 pray(median, vec);
    Ray2 nray(median, -vec);

    std::vector<Point2> intersections{};
    /* iterate over polygon edges to find intersections with separating line */
    for (auto eit = poly.edges_begin(); eit != poly.edges_end(); ++eit) {
        auto s = *eit;
        // Compute intersection with *positive* ray
        auto inter = CGAL::intersection(s, pray);

        // if we find an intersection, then add the intersection point to the list of intersections
        if (inter.has_value())
            intersections.push_back(boost::get<Point2>(inter.value()));
        inter = CGAL::intersection(s, nray);
        if (inter.has_value())
            intersections.push_back(boost::get<Point2>(inter.value()));
    }

    // remove duplicates
    auto last = distinct(intersections.begin(), intersections.end(), P2Comp{});
    intersections.erase(last, intersections.end());

    std::vector<Point2> b1(intersections.begin(), intersections.end()), b2(intersections.begin(), intersections.end());

    // add to left or right
    for (auto eit = poly.edges_begin(); eit != poly.edges_end(); ++eit) {
        Segment2 s = *eit;
        add_to_bisection(b1, b2, vec, median, s.source());
        add_to_bisection(b1, b2, vec, median, s.target());
    }

    std::vector<Point2> chull1{}, chull2{};

    CGAL::ch_graham_andrew(b1.begin(), b1.end(), std::back_inserter(chull1));
    CGAL::ch_graham_andrew(b2.begin(), b2.end(), std::back_inserter(chull2));

    Polygon2 poly1(chull1.begin(), chull1.end());
    Polygon2 poly2(chull2.begin(), chull2.end());

    return {poly1, poly2};
}
