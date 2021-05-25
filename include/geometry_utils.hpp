//
// Created by xetql on 3/17/21.
//

#pragma once

#include <vector>
#include <array>
#include <cassert>
#include "serial/format.hpp"
#include "parallel/algorithm.hpp"
#include <ostream>
#include "algorithm.hpp"
#include "numeric/utils.hpp"
#include <sstream>
#include <mpi.h>
namespace norcb{
const double acceptable_error = 1e-12;

template<class Real, class ForwardIt, class GetVelocityFunc>
std::pair<Real, Real> par_get_average_velocity(ForwardIt el_begin, ForwardIt el_end, unsigned longest_axis, MPI_Comm comm,
                                               GetVelocityFunc getVelocity) {

    auto size = std::distance(el_begin, el_end);

    std::array<Real, 2> s = {0., 0.};
    auto folding_axis = ((longest_axis + 1) % 2);

    decltype(size) total_size;

    MPI_Allreduce(&size, &total_size, 1, par::get_mpi_type<decltype(size)>(), MPI_SUM, comm);

    // O(n/p)
    for (auto i = 0; i < size; ++i) {
        auto velocity = *getVelocity(&(*(el_begin + i)));
        const auto d2 = velocity[0]*velocity[0] + velocity[1]*velocity[1];
        velocity[0] /= d2; velocity[1] /= d2;
        if(velocity.at(folding_axis) < 0) {
            s[0] += -velocity.at(0);
            s[1] += -velocity.at(1);
        } else {
            s[0] += velocity.at(0);
            s[1] += velocity.at(1);
        }
    }

    MPI_Allreduce(MPI_IN_PLACE, s.data(), 2, par::get_mpi_type<Real>() , MPI_SUM, comm);

    return {s[0] / total_size, s[1] / total_size};
}



template<class Real>
Real side(Real tx, Real ty, Real dx, Real dy) {
    return dx * ty - dy * tx;
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



inline std::pair<std::vector<double>, std::vector<double>> rotate_copy(const std::array<std::array<double, 2>, 2> &matrix, const std::vector<double> &x,
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
}

