//
// Created by xetql on 3/17/21.
//

#pragma once

#include <vector>
#include <array>
#include <cassert>
#include <serial/format.hpp>
#include <ostream>
#include <parallel/algorithm.hpp>
#include "algorithm.hpp"
#include "numeric/utils.hpp"
#include <sstream>

static const double acceptable_error = 1e-12;

struct Point {
    long double x, y;

    Point operator-(const Point& p) const {
        return Point {x-p.x, y-p.y};
    }
    Point operator+(const Point& p) const {
        return Point {x+p.x, y+p.y};
    }
    Point operator/(const Point& p) const {
        return Point {x/p.x, y/p.y};
    }
    Point operator*(const Point& p) const {
        return Point {x*p.x, y*p.y};
    }
    Point operator*(const long double& v) const {
        return Point {x*v, y*v};
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

    Segment(const Point& a, const Point& b) : a(a), b(b){}

    Segment operator-(const Point &p) const {
        return Segment { Point {a.x - p.x, a.y - p.y}, Point {b.x - p.x, b.y - p.y} };
    }

    [[nodiscard]] bool contains (const Point& p) const {

        const auto crossproduct = (p.y - a.y) * (b.x - a.x) - (p.x - a.x) * (b.y - a.y);

        if(std::fabs(crossproduct) > acceptable_error) {
            par::pcout() << *this << " doesnt contains " << p << " due to crossproduct: " << crossproduct << std::endl;
            return false;
        }

        const auto dotproduct = (p.x - a.x) * (b.x - a.x) + (p.y - a.y) * (b.y - a.y);

        if(dotproduct < acceptable_error) {
            par::pcout() << *this << " doesnt contains " << p << " due to dotproduct: " << dotproduct << std::endl;

            return false;
        }

        const auto squaredlengthab = (b.x - a.x) * (b.x - a.x) + (b.y - a.y) * (b.y - a.y);

        if(dotproduct > squaredlengthab) {
            par::pcout() << *this << " doesnt contains " << p << " due to dotproduct: " << dotproduct << std::endl;
            return false;
        }

        return true;
    }

    friend std::ostream &operator<<(std::ostream &os, const Segment &segment) {
        os << "a: " << segment.a << " b: " << segment.b;
        return os;
    }
};

inline auto sqr_dist( const Segment& s, const Point& p )
{
    auto  n = s.b - s.a;
    auto pa = s.a - p;
    auto bp = p - s.b;

    auto c = n.x*pa.x + n.y*pa.y; //Dot( n, pa );

    // Closest point is a
    if ( c > 0.0 )
        return pa.x*pa.x + pa.y*pa.y; // Dot( pa, pa );


    // Closest point is b
    c = n.x*bp.x + n.y*bp.y;
    if ( c > 0.0 )
        return bp.x*bp.x + bp.y*bp.y;

    // Closest point is between a and b
    auto e = pa - n * (c / (n.x*n.x + n.y*n.y ));

    return e.x*e.x + e.y*e.y;
}

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

inline auto crossproduct(const Point& p, const Point& q, const Point& r) {
    return (q.x * r.y - r.x * q.y) - (p.x * r.y - r.x * p.y) + (p.x * q.y - p.y * q.x);
}

inline int orientation_sign(const Point& p, const Point& q, const Point& r)
{
    const auto res = crossproduct(p, q, r);
    return (almost_equal(res, 0.0l, 4)) ? 0 : ((res > 0.0) ? 1 : -1);
}

class Polygon {
    template<class FIt> void create_polygon(FIt beg, FIt end) {
        auto n = std::distance(beg, end);

        std::vector<Point> pUpper{}; pUpper.reserve(n);
        std::vector<Point> pLower{}; pLower.reserve(n);

        //sort by x axis
        std::sort(beg, end, [](auto& p1, auto& p2) {
            return ((almost_equal(p1.x, p2.x, 4)) ? (p1.y < p2.y) : (p1.x < p2.x));
        });

        auto& right = beg[n - 1];
        auto& left  = beg[0];

        pLower.push_back(left);

        for (auto i = 1; i < n - 1; ++i) {
            if (orientation_sign(right, left, beg[i]) > 0) {
                pLower.push_back(beg[i]);
            } else {
                pUpper.push_back(beg[i]);
            }
        }
        pUpper.push_back(right);
        std::reverse(pUpper.begin(), pUpper.end());

        vertices.assign(pLower.begin(), pLower.end());
        vertices.insert(vertices.end(), pUpper.begin(), pUpper.end());

        edges.reserve(n);
        for(auto i = 1; i <= n; ++i){
            const auto& p1 = vertices.at(i - 1);
            const auto& p2 = vertices.at(i % n);
            edges.emplace_back(p1, p2);
        }
    }
public:
    std::vector<Point> vertices{};
    std::vector<Segment> edges{};
    Polygon(){}

    template<class FIt> Polygon(FIt beg, FIt end){
        create_polygon(beg, end);
    }

    [[nodiscard]] auto edges_cbegin() const {
        return edges.cbegin();
    }
    [[nodiscard]] auto edges_cend() const {
        return edges.cend();
    }
    [[nodiscard]] auto edges_begin() const {
        return edges.begin();
    }
    [[nodiscard]] auto edges_end() const {
        return edges.end();
    }
    [[nodiscard]] unsigned size() const {
        return vertices.size();
    }
};

inline std::optional<Point> intersection(const Line& l1, const Line& l2) {
    const auto&[p1, p2] = l1.get_points();
    const auto&[p3, p4] = l2.get_points();

    const auto D = (p1.x-p2.x)*(p3.y-p4.y) - (p1.y-p2.y)*(p3.x-p4.x);

    if(std::fabs(D) < std::numeric_limits<double>::epsilon()) {
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
    const auto s = sign(side(v.x, v.y, median_to_p.x, median_to_p.y));

    if (s <= 0) {
        b1.push_back(p);
    } else {
        b2.push_back(p);
    }
}

inline auto get_angle(const Vector &v, const Vector &origin) {
    auto dot = v.x * origin.x + v.y * v.y;
    auto det = v.x * origin.y - v.y * origin.x;
    return std::atan2(det, dot);
}

template<class Real>
inline auto get_angle(Real x1, Real y1, Real x2, Real y2) {
    auto dot = x1 * x2 + y1 * y2;
    auto det = x1 * y2 - y1 * x2;
    return std::atan2(det, dot);
}



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
/*
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
*/
template<class Real>
std::pair<Polygon, Polygon> bisect_polygon(const Polygon &poly, Real vx, Real vy, Real px, Real py) {

    Vector vec {vx, vy};
    Point  pmed{px, py};
    Line   med(pmed, vec);

    std::vector<Point> intersections {};

    /* iterate over polygon edges to find intersections with separating line...*/

    for (auto eit = poly.edges_cbegin(); eit != poly.edges_cend(); eit++) {
        const auto& s = *eit;
        // Compute intersection
        const auto opt = intersection(med, s);

        // if we find an intersection, then add the intersection point to the list of intersections
        if (opt.has_value()) {
            const auto& inter = opt.value();
            intersections.push_back(inter);
        }
    }

    if (intersections.size() != 2) {
        std::stringstream ss;
        std::for_each(intersections.begin(), intersections.end(), [&ss](auto p) {ss << p;});
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if(!rank) throw std::logic_error(ser::fmt("A line intersects 2 polygon edges. %s", ss.str()));
    }

    std::vector<Point> b1(intersections.begin(), intersections.end()), b2(intersections.begin(), intersections.end());
    b1.reserve(b1.size() + poly.size());
    b2.reserve(b2.size() + poly.size());

    // add to left or right
    for (auto eit = poly.edges_cbegin(); eit != poly.edges_cend(); eit++) {
        const auto& s = *eit;
        add_to_bisection(b1, b2, vec, pmed, s.a);
        add_to_bisection(b1, b2, vec, pmed, s.b);
    }

    Polygon poly1(b1.begin(), b1.end());
    Polygon poly2(b2.begin(), b2.end());

    return {poly1, poly2};
}
