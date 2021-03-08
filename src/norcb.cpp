//
// Created by xetql on 12/27/20.
//

#include "norcb.hpp"
#include <CGAL/iterator.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/ch_graham_andrew.h>
#include <boost/optional/optional_io.hpp>

std::array<std::array<double, 2>, 2> get_rotation_matrix(double theta) {
    std::array<std::array<double, 2>, 2> rotate{};
    rotate[0][0] = std::cos(theta);
    rotate[0][1] = -std::sin(theta);
    rotate[1][0] = std::sin(theta);
    rotate[1][1] = std::cos(theta);
    return rotate;
}

Polygon2 init_domain(Real minx, Real miny, Real maxx, Real maxy){
    Point2 p1(minx, miny), p2(maxx, miny), p3(maxx, maxy), p4(minx, maxy);
    Polygon2 d;
    d.push_back(p1);
    d.push_back(p2);
    d.push_back(p3);
    d.push_back(p4);
    return d;
}

template<class T>
typename std::enable_if<!std::numeric_limits<T>::is_integer, bool>::type almost_equal(T x, T y, int ulp) {
    // the machine epsilon has to be scaled to the magnitude of the values used
    // and multiplied by the desired precision in ULPs (units in the last place)
    return std::fabs(x-y) <= std::numeric_limits<T>::epsilon() * std::fabs(x+y) * ulp
           // unless the result is subnormal
           || std::fabs(x-y) < std::numeric_limits<T>::min();
}
/**
 * Distinct without sorting worst is O(N^2) and best is O(N)
 * @tparam ForwardIt
 * @tparam BinaryComp
 * @param first
 * @param last
 * @param eq
 * @return
 */
template<class ForwardIt, class BinaryComp>
ForwardIt distinct(ForwardIt first, ForwardIt last, BinaryComp eq) {
    if(first == last) return last;
    while(first != last) {
        ForwardIt next = std::next(first);
        while(next != last) {
            if(eq(*first, *next)) { // not alone, swap and erase
                *next = *(last-1);
                last--;
            } else next++;
        }
        first++;
    }
    return last;
}

struct P2Comp{
    int operator()(Point2 p1, Point2 p2) {
        return almost_equal(p1.x(), p2.x(), 2) && almost_equal(p1.y(), p2.y(), 2);
    }
};

double side(const Vector2& T, const Vector2& D){
    return D.x() * T.y() - D.y() * T.x();
}

std::vector<int> side(const Vector2& T, const Point2& p, const std::vector<double>& x, const std::vector<double>& y){
    const auto size = x.size();
    std::vector<int> sides(size, 0);
    const double px = p.x(), py = p.y();
    for(unsigned i = 0; i < size; ++i) {
        sides[i] = (x[i] - px) * T.y() - (y[i] - py) *T.x();
    }
    return sides;
}

int sign(double x) {
    return x < 0 ? -1 : (x == 0 ? 0 : 1);
}

// def get_angle(v, origin):
// dot = np.dot(v, origin)
// det = np.linalg.det([v, origin])
// return np.arctan2(det, dot)


// def get_average_velocity(velocities, axis):
// v = np.copy(velocities)
// R = get_rotation_matrix(np.pi)
// # print(np.inner(R, v))
// msk = np.where(v[:, axis] >= 0)
// v[msk] = np.inner(R, v[msk]).T
// return np.mean(v, axis=0)

Vector2 compute_average_velocity(std::vector<double> vx, std::vector<double> vy, bool axis){
    const auto size = vx.size();
    auto R = get_rotation_matrix(M_PI);
    std::vector<double>& target = !axis ? vx : vy;
    for(unsigned i = 0; i < size; ++i){
        target[i] = target[i] >= 0 ? -target[i] : target[i];
    }
    return Vector2(std::accumulate(vx.begin(), vx.end(), 0.)/size,
                   std::accumulate(vy.begin(), vy.end(), 0.)/size);
}

double get_angle(const Vector2& v, const Vector2& origin){
    auto dot = v * origin;
    auto det = v.x()*origin.y() - v.y()*origin.x();
    return std::atan2(det, dot);
}

std::pair<std::vector<double>, std::vector<double>>
rotate_copy(const std::array<std::array<double, 2>, 2>& matrix, const std::vector<double>& x, const std::vector<double>& y){
    const auto size = x.size();
    std::vector<double> rx(size, 0), ry(size, 0);
#pragma GCC ivdep
    for(unsigned i = 0; i < size; ++i) {
        rx[i] = matrix[0][0]*x[i] + matrix[0][1]*y[i];
        ry[i] = matrix[1][0]*x[i] + matrix[1][1]*y[i];
    }
    return {rx, ry};
}

std::pair<double, double> rotate(const std::array<std::array<double, 2>, 2>& matrix, double x, double y){
    return {matrix[0][0]*x + matrix[0][1]*y, matrix[1][0]*x + matrix[1][1]*y};
}

void rotate(const std::array<std::array<double, 2>, 2>& matrix, std::vector<double>& x, std::vector<double>& y){
    const auto size = x.size();
    std::vector<double> rx(size, 0), ry(size, 0);
    double rxi, ryi;
#pragma GCC ivdep
    for(unsigned i = 0; i < size; ++i) {
        rxi = matrix[0][0]*x[i] + matrix[0][1]*y[i];
        ryi = matrix[1][0]*x[i] + matrix[1][1]*y[i];
        x[i] = rxi;
        y[i] = ryi;
    }

}

void add_to_bisection(std::vector<Point2>& b1, std::vector<Point2>& b2, const Vector2& v, const Point2& pmed, const Point2& p){
    auto s = sign(side(v, p - pmed));
    if (s <= 0) {
        b1.push_back(p);
    } else {
        b2.push_back(p);
    }
}

std::pair<Polygon2, Polygon2> bisect_polygon(const Polygon2& poly, const Vector2& vec, const Point2& median) {
    Ray2 pray(median,  vec);
    Ray2 nray(median, -vec);

    std::vector<Point2> intersections{};
    for(auto eit = poly.edges_begin(); eit != poly.edges_end(); ++eit){
        auto s = *eit;
        auto inter = CGAL::intersection(s, pray);
        if (inter.has_value())
            intersections.push_back(boost::get<Point2>(inter.value()));
        inter = CGAL::intersection(s, nray);
        if (inter.has_value())
            intersections.push_back(boost::get<Point2>(inter.value()));
    }

    auto last = distinct(intersections.begin(), intersections.end(), P2Comp {});
    intersections.erase(last, intersections.end());

    std::vector<Point2> b1(intersections.begin(), intersections.end()), b2(intersections.begin(), intersections.end());
    for(auto eit = poly.edges_begin(); eit != poly.edges_end(); ++eit) {
        Segment2 s = *eit;
        add_to_bisection(b1, b2, vec, median, s.source());
        add_to_bisection(b1, b2, vec, median, s.target());
    }

    std::vector<Point2> chull1, chull2;

    CGAL::ch_graham_andrew(b1.begin(), b1.end(), std::back_inserter(chull1));
    CGAL::ch_graham_andrew(b2.begin(), b2.end(), std::back_inserter(chull2));

    Polygon2 poly1(chull1.begin(), chull1.end());
    Polygon2 poly2(chull2.begin(), chull2.end());

    return {poly1, poly2};
}

NoRCB partition(Integer P, std::vector<double>& x,  std::vector<double>& y,
                           std::vector<double>& vx, const std::vector<double>& vy,
                           const Polygon2& domain) {
    const auto size = x.size();
    std::vector<std::tuple<Polygon2,
        std::vector<double>, std::vector<double>,
        std::vector<double>, std::vector<double>>> partitions;
    partitions.emplace_back(domain, x, y, vx, vy);
    while(P > 1) {
        decltype(partitions) bisected_parts{};
        for(auto& partition : partitions){
            auto& [domain, x, y, vx, vy] = partition;
            auto [minx, maxx] = std::minmax(x.begin(), x.end());
            auto [miny, maxy] = std::minmax(x.begin(), x.end());

            unsigned target_axis;

            if((maxx-minx) > (maxy-miny)){
                target_axis = 0;
            } else {
                target_axis = 1;
            }

            Vector2 origin(0., 1.);
            auto avg_vel = compute_average_velocity(vx, vy, target_axis);

            auto theta         = get_angle(avg_vel, origin);
            auto clockwise     = get_rotation_matrix(theta);
            auto anticlockwise = get_rotation_matrix(-theta);

            rotate(clockwise, x, y);

            avg_vel = Vector2(
                clockwise[0][0] * avg_vel.x() + clockwise[0][1]*avg_vel.y(),
                clockwise[1][0] * avg_vel.x() + clockwise[1][1]*avg_vel.y()
            );

            double norm = std::sqrt(avg_vel.x()*avg_vel.x() + avg_vel.y()*avg_vel.y());
            avg_vel = avg_vel / norm;

            double median;
            {
                std::vector<double> rx(x.begin(), x.end());
                std::nth_element(rx.begin(), rx.begin() + rx.size() / 2, rx.end());
                median = *(rx.begin() + rx.size() / 2);
            }

            auto c = std::count_if(x.begin(), x.end(), [median](auto v){return v <= median;});

            std::vector<double> lx(c), ly(c), rx(size-c), ry(size-c);
            std::vector<double> lvx(c), lvy(c), rvx(size-c), rvy(size-c);
            unsigned l = 0, r = 0;
            for(unsigned i = 0; i < size; ++i) {
                if(x[i] <= median) {
                    lx[l]  = x[i];
                    ly[l]  = y[i];
                    lvx[l] = vx[i];
                    lvy[l] = vy[i];
                    l++;
                } else {
                    rx[l]  = x[i];
                    ry[l]  = y[i];
                    rvx[l] = vx[i];
                    rvy[l] = vy[i];
                    r++;
                }
            }
            auto [pmedx, pmedy] = rotate(anticlockwise, median, 1.0);
            Point2 pmedian(pmedx, pmedy);

            rotate(anticlockwise, lx, ly);
            rotate(anticlockwise, rx, ry);
            auto [lpoly, rpoly] = bisect_polygon(domain, avg_vel, pmedian);
            bisected_parts.emplace_back(lpoly, lx, ly, lvx, lvy);
            bisected_parts.emplace_back(rpoly, rx, ry, rvx, rvy);
        }
        P /= 2;
        partitions = bisected_parts;
    }

    return NoRCB();
}
