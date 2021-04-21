//
// Created by xetql on 3/17/21.
//
#include "geom/utils.hpp"

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
    auto s = sign(CGAL::to_double(side(v, p - pmed)));
    if (s <= 0) {
        b1.push_back(p);
    } else {
        b2.push_back(p);
    }
}

