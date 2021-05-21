//
// Created by xetql on 3/17/21.
//
#include "geom/utils.hpp"


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



