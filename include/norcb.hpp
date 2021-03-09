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

struct NoRCB {};

Polygon2 init_domain(Real minx, Real miny, Real maxx, Real maxy);

std::pair<Polygon2, Polygon2> bisect_polygon(const Polygon2& poly, const Vector2& vec, const Point2& median);

NoRCB partition(Integer P, const std::vector<Point2>&, const Polygon2& domain);

#endif //NORCB_NORCB_HPP
