//
// Created by xetql on 12/27/20.
//

#ifndef NORCB_NORCB_HPP
#define NORCB_NORCB_HPP

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Polygon_2.h>
#include <CGAL/Point_2.h>
#include <CGAL/Vector_2.h>
#include <CGAL/Ray_2.h>

#include <tuple>
#include <iostream>

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Polygon2 = CGAL::Polygon_2<K>;
using Segment2 = CGAL::Segment_2<K>;
using Point2 = CGAL::Point_2<K>;
using Vector2 = CGAL::Vector_2<K>;
using Ray2 = CGAL::Ray_2<K>;

using Integer = unsigned;
using Real    = float;

namespace norcb {

struct NoRCB {
    friend std::ostream &operator<<(std::ostream &os, const NoRCB &rcb);
};
Polygon2 init_domain(Real minx, Real miny, Real maxx, Real maxy);

std::pair<Polygon2, Polygon2> bisect_polygon(const Polygon2& poly, const Vector2& vec, const Point2& median);
namespace seq {
std::vector<std::tuple<Polygon2,
std::vector<double>, std::vector<double>,
std::vector<double>, std::vector<double>>> partition(Integer P, std::vector<double>& x, std::vector<double>& y,
                                                     std::vector<double>& vx, std::vector<double>& vy,
                                                             const Polygon2& domain);
}
namespace par {
std::tuple<Polygon2, std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>>
partition(Integer P, const std::vector<double>& x, const std::vector<double>& y,
          const std::vector<double>& vx, const std::vector<double>& vy,
          const Polygon2& domain);
}
}
#endif //NORCB_NORCB_HPP
