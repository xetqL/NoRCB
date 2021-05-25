//
// Created by xetql on 5/25/21.
//

#ifndef YALBB_GEOMETRY_HPP
#define YALBB_GEOMETRY_HPP

#include <utility>
#include <geometry_utils.hpp>
#include "numeric/utils.hpp"

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
            return false;
        }

        const auto dotproduct = (p.x - a.x) * (b.x - a.x) + (p.y - a.y) * (b.y - a.y);

        if(dotproduct < -acceptable_error) {
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
inline int orientation_sign(const Point& p, const Point& q, const Point& r) {
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
        if(s.contains(i)){
            auto& p = i;
            auto[a,b] = s;
            auto pa = a-p;
            auto pb = b-p;
            auto lpa = pa.x*pa.x+pa.y*pa.y;
            auto lpb = pb.x*pb.x+pb.y*pb.y;
            if(lpa < acceptable_error || lpb < acceptable_error) {
                if(lpa < lpb) return a;
                else return b;
            }
            return i;
        } else {
            return std::nullopt;
        }

    } else {
        return std::nullopt;
    }
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
Polygon init_domain(Real minx, Real miny, Real maxx, Real maxy) {
    std::array<Point, 4> p = { Point{minx, miny},Point{maxx, miny},Point{maxx, maxy},Point{minx, maxy} };
    Polygon d(p.begin(), p.end());
    return d;
}
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
    auto last = distinct(intersections.begin(), intersections.end(), [](auto& a, auto& b){
        auto ba = a-b;
        return std::sqrt(ba.x*ba.x+ba.y*ba.y) < 1e-12;
    });
    intersections.erase(last, intersections.end());

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

#endif //YALBB_GEOMETRY_HPP
