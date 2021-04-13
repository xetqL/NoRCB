//
// Created by xetql on 12/27/20.
//

#ifndef NORCB_NORCB_HPP
#define NORCB_NORCB_HPP

#include <parallel/algorithm.hpp>

#include "geom/utils.hpp"

#include <boost/optional/optional_io.hpp>

#include <tuple>
#include <iostream>

#include <mpi.h>

using Integer = unsigned;
using Real    = float;

template<class T> std::ostream &operator<<(std::ostream &os, const std::vector<T> &v) {
    std::for_each(v.begin(), v.end() - 1, [&os](auto &x) { os << x << ","; });
    os << *(v.end() - 1) << std::endl;
    return os;
}

namespace norcb {
namespace{
template<class ReduceStrategy>
Vector2 compute_average_velocity(std::vector<double> vx, std::vector<double> vy, bool axis, ReduceStrategy reduce) {
    const auto size = vx.size();
    auto R = get_rotation_matrix(M_PI);
    std::vector<double> &target = !axis ? vx : vy;

    for (unsigned i = 0; i < size; ++i) {
        const auto rvx = R[0][0] * vx[i] + R[0][1] * vy[i];
        const auto rvy = R[1][0] * vx[i] + R[1][1] * vy[i];
        vx[i] = target[i] >= 0 ? rvx : vx[i];
        vy[i] = target[i] >= 0 ? rvy : vy[i];
    }
    const auto[acc_vx, xsize] = reduce(vx.begin(), vx.end());
    const auto[acc_vy, ysize] = reduce(vy.begin(), vy.end());


    return Vector2(acc_vx / xsize, acc_vy / ysize);
}
}

struct NoRCB {
    K trait {};
    int world_size, rank;

    std::vector<Polygon2> subdomains {};
    MPI_Comm comm;

    NoRCB(std::vector<Polygon2> subdomains, MPI_Comm comm) : comm(comm), subdomains(std::move(subdomains)) {
        MPI_Comm_size(comm, &world_size);
        MPI_Comm_rank(comm, &rank);
    }

    template<class Real>
    bool is_inside(Polygon2& domain, Real x, Real y) {
        Point2 pt(x,y);
        return CGAL::bounded_side_2(domain.vertices_begin(), domain.vertices_end(), pt, trait) != CGAL::ON_UNBOUNDED_SIDE;
    }

    template<class Real>
    void get_owner(Real x, Real y, int* owner){
        for((*owner) = 0; (*owner) < world_size; (*owner)++){
            if(is_inside(subdomains.at(*owner), x, y)) {
                return;
            }
        }
        *owner = -1;
    }

    template<class Real>
    void get_intersecting_domains(Real x1, Real x2, Real y1, Real y2, Real z1, Real z2, int* PEs, int* num_found) {
        Point2 x1y1 = Point2(x1, y1),
               x1y2 = Point2(x1, y2),
               x2y2 = Point2(x2, x2),
               x2y1 = Point2(x2, y1);

        Polygon2 interaction_square;
        interaction_square.push_back(x1y1);
        interaction_square.push_back(x1y2);
        interaction_square.push_back(x2y2);
        interaction_square.push_back(x2y1);

        *num_found = 0;
        for(auto PE = 0; PE < world_size; ++PE){
            const Polygon2& subdomain = subdomains.at(PE);
            if(CGAL::do_intersect(interaction_square, subdomain)){
                PEs[*num_found] = PE;
                *num_found += 1;
            }
        }
    }

};

NoRCB* allocate_from(NoRCB* from);
void destroy(NoRCB* lb);

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

namespace seq {
namespace{
std::vector<std::tuple<Polygon2,
std::vector<double>, std::vector<double>,
std::vector<double>, std::vector<double>>> partition(Integer P, std::vector<double>& x, std::vector<double>& y,
                                                     std::vector<double>& vx, std::vector<double>& vy,
                                                             const Polygon2& domain);
}
}

namespace parallel {
namespace {
Vector2 compute_average_velocity(std::vector<double> vx, std::vector<double> vy, bool axis, MPI_Comm comm);

template<class ForwardIt>
Vector2
parallel_compute_average_velocity(ForwardIt vx_begin, ForwardIt vx_end, ForwardIt vy_begin, bool axis, MPI_Comm comm) {
    auto size = std::distance(vx_begin, vx_end);
    auto R = get_rotation_matrix(M_PI);
    ForwardIt target_begin = !axis ? vx_begin : vy_begin;
    std::array<double, 2> s = {0., 0.};

    for (unsigned i = 0; i < size; ++i) {
        const auto rvx = R[0][0] * (*(vx_begin + i)) + R[0][1] * (*(vy_begin + i));
        const auto rvy = R[1][0] * (*(vx_begin + i)) + R[1][1] * (*(vy_begin + i));
        s[0] += (*(target_begin + i)) >= 0 ? rvx : (*(vx_begin + i));
        s[1] += (*(target_begin + i)) >= 0 ? rvy : (*(vy_begin + i));
    }

    MPI_Allreduce(MPI_IN_PLACE, s.data(), 2, MPI_DOUBLE, MPI_SUM, comm);
    MPI_Allreduce(MPI_IN_PLACE, &size, 1, par::get_mpi_type<decltype(size)>(), MPI_SUM, comm);

    return Vector2(s[0] / size, s[1] / size);
}

template<class Real, class ForwardIt, class GetPositionFunc, class GetVelocityFunc>
Vector2 parallel_compute_average_velocity(ForwardIt el_begin, ForwardIt el_end, unsigned target_axis, MPI_Comm comm,
                                          GetPositionFunc getPosition, GetVelocityFunc getVelocity) {
    auto size = std::distance(el_begin, el_end);
    auto R = get_rotation_matrix(M_PI);
    std::array<Real, 2> s = {0., 0.};

    for (unsigned i = 0; i < size; ++i) {
        auto position = getPosition(&(*(el_begin + i)));
        auto velocity = getVelocity(&(*(el_begin + i)));

        const auto rvx = R[0][0] * (velocity->at(0)) + R[0][1] * (velocity->at(1));
        const auto rvy = R[1][0] * (velocity->at(0)) + R[1][1] * (velocity->at(1));

        s[0] += (position->at(target_axis)) >= 0 ? rvx : (velocity->at(0));
        s[1] += (position->at(target_axis)) >= 0 ? rvy : (velocity->at(1));
    }

    MPI_Allreduce(MPI_IN_PLACE, s.data(), 2, par::get_mpi_type<Real>() , MPI_SUM, comm);
    MPI_Allreduce(MPI_IN_PLACE, &size, 1, par::get_mpi_type<decltype(size)>(), MPI_SUM, comm);

    return Vector2(s[0] / size, s[1] / size);
}
}
std::vector<std::tuple<Polygon2,
        std::vector<double>, std::vector<double>,
        std::vector<double>, std::vector<double>>> partition(Integer P, std::vector<double>& x, std::vector<double>& y,
                                                             std::vector<double>& vx, std::vector<double>& vy,
                                                             const Polygon2& domain, MPI_Comm comm);

template<class ForwardIt>
std::vector<std::tuple<Polygon2, ForwardIt, ForwardIt, ForwardIt, ForwardIt, ForwardIt>>
partition(Integer P,
          ForwardIt x_begin, ForwardIt x_end,
          ForwardIt y_begin, ForwardIt vx_begin, ForwardIt vy_begin,
          const Polygon2 &domain, MPI_Comm comm) {

    std::vector<std::tuple<Polygon2, ForwardIt, ForwardIt, ForwardIt, ForwardIt, ForwardIt>> partitions {};

    partitions.emplace_back(domain, x_begin, x_end, y_begin, vx_begin, vy_begin);

    unsigned npart = partitions.size();
    while (npart != P) {
        decltype(partitions) bisected_parts{};
        for (auto &partition : partitions) {
            auto&[domain, x_begin, x_end, y_begin, vx_begin, vy_begin] = partition;
            const unsigned elements_in_subdomain = std::distance(x_begin, x_end);
            const ForwardIt y_end  = y_begin  + elements_in_subdomain,
                      vx_end = vx_begin + elements_in_subdomain,
                      vy_end = vy_begin + elements_in_subdomain;

            const double minx = *std::min_element(x_begin, x_end);
            const double maxx = *std::max_element(x_begin, x_end);
            const double miny = *std::min_element(y_begin, y_end);
            const double maxy = *std::max_element(y_begin, y_end);

            const unsigned target_axis = ((maxx - minx) > (maxy - miny)) ? 0 : 1;

            const Vector2 origin(0., 1.);
            auto avg_vel = parallel_compute_average_velocity(vx_begin, vx_end, vy_begin, target_axis, comm);

            const auto theta         = get_angle(avg_vel, origin);
            const auto clockwise     = get_rotation_matrix(theta);
            const auto anticlockwise = get_rotation_matrix(-theta);

            rotate(clockwise, x_begin, x_end, y_begin);

            const double norm = std::sqrt(avg_vel.x() * avg_vel.x() + avg_vel.y() * avg_vel.y());
            avg_vel = avg_vel / norm;

            double median;
            {
                std::remove_const_t<decltype(elements_in_subdomain)>  size {0};
                MPI_Allreduce(&elements_in_subdomain, &size, 1, par::get_mpi_type<decltype(size)>(), MPI_SUM, comm);
                auto midIdx = size / 2;
                median = par::find_nth(x_begin, x_end, midIdx, comm);
            }

            const auto c = std::count_if(x_begin, x_end, [median](auto v) { return v <= median; });

            ForwardIt x_median  = x_begin  + c;
            ForwardIt y_median  = y_begin  + c;
            ForwardIt vx_median = vx_begin + c;
            ForwardIt vy_median = vy_begin + c;

            unsigned l = 0, r = 0, i = 0;
            while(l < c) {
                if (*(x_begin + i) <= median) {
                    l++;
                }
                if (*(x_begin + i) > median) {
                    std::iter_swap((x_begin +  l),  (x_begin + c + r));
                    std::iter_swap((y_begin +  l),  (y_begin + c + r));
                    std::iter_swap((vx_begin + l), (vx_begin + c + r));
                    std::iter_swap((vy_begin + l), (vy_begin + c + r));
                    r++;
                }
            }

            rotate(anticlockwise, x_begin, x_median, y_begin);

            const auto[pmedx, pmedy] = rotate(anticlockwise, median, 1.0);

            const Point2 pmedian(pmedx, pmedy);

            rotate(anticlockwise, x_median, x_end, y_median);

            const auto[lpoly, rpoly] = bisect_polygon(domain, avg_vel, pmedian);

            bisected_parts.emplace_back(lpoly, x_begin,  x_median, y_begin,  vx_begin,  vy_begin);
            bisected_parts.emplace_back(rpoly, x_median, x_end,    y_median, vx_median, vy_median);
        }
        npart *= 2;
        partitions = bisected_parts;
    }

    return partitions;
}


template<class Real, class ForwardIt, class GetPositionFunc, class GetVelocityFunc>
std::vector<std::tuple<Polygon2, ForwardIt, ForwardIt>>
partition(Integer P, ForwardIt el_begin, ForwardIt el_end,
          const Polygon2 &domain, MPI_Datatype datatype, MPI_Comm comm,
          GetPositionFunc getPosition, GetVelocityFunc getVelocity) {

    std::vector<std::tuple<Polygon2, ForwardIt, ForwardIt>> partitions {};

    partitions.emplace_back(domain, el_begin, el_end);

    unsigned npart = partitions.size();
    while (npart != P) {
        decltype(partitions) bisected_parts{};
        for (auto &partition : partitions) {
            auto&[domain, el_begin, el_end] = partition;
            const unsigned elements_in_subdomain = std::distance(el_begin, el_end);

            Real minx = std::numeric_limits<Real>::max(),
                maxx = std::numeric_limits<Real>::lowest(),
                miny = std::numeric_limits<Real>::max(),
                maxy = std::numeric_limits<Real>::lowest();

            std::for_each(el_begin, el_end, [&minx, &miny, &maxx, &maxy, getPosition] (const auto& e1) {
                auto position = getPosition(&e1);
                minx = std::min(position->at(0), minx);
                maxx = std::max(position->at(0), maxx);
                miny = std::min(position->at(1), miny);
                maxy = std::max(position->at(1), maxy);
            });

            const unsigned target_axis = ((maxx - minx) > (maxy - miny)) ? 0 : 1;

            const Vector2 origin(0., 1.);
            auto avg_vel = parallel_compute_average_velocity<Real>(el_begin, el_end, target_axis, comm, getPosition, getVelocity);

            const auto theta         = get_angle(avg_vel, origin);
            const auto clockwise     = get_rotation_matrix(theta);
            const auto anticlockwise = get_rotation_matrix(-theta);

            rotate(clockwise, el_begin, el_end, getPosition);

            const auto norm = std::sqrt(avg_vel.x() * avg_vel.x() + avg_vel.y() * avg_vel.y());
            avg_vel = avg_vel / norm;

            Real median;
            {
                std::remove_const_t<decltype(elements_in_subdomain)>  size {0};
                MPI_Allreduce(&elements_in_subdomain, &size, 1, par::get_mpi_type<decltype(size)>(), MPI_SUM, comm);
                auto midIdx = size / 2-1;
                auto el_median = par::find_nth(el_begin, el_end, midIdx, datatype, comm, [getPosition](const auto& a, const auto& b){
                    auto posa = getPosition(&a);
                    auto posb = getPosition(&b);
                    return posa->at(0) < posb->at(0);
                }, [getPosition](const auto& a, const auto& b){
                    auto posa = getPosition(&a);
                    auto posb = getPosition(&b);
                    return posa->at(0) == posb->at(0);
                });
                median = getPosition(&el_median)->at(0);
            }

            const auto c = std::count_if(el_begin, el_end, [getPosition, median](const auto& v) {
                auto pos = getPosition(&v);
                return pos->at(0) <= median;
            });
            ForwardIt el_median = el_begin + c;

            unsigned l = 0, r = 0, i = 0;
            while(l < c) {
                auto position = getPosition(&(*(el_begin + i)));
                if (position->at(0) <= median) {
                    l++;
                }
                if (position->at(0) > median) {
                    std::iter_swap(el_begin + i, el_begin + c + r);
                    r++;
                }
            }

            rotate(anticlockwise, el_begin, el_end, getPosition);

            const auto[pmedx, pmedy] = rotate(anticlockwise, median, 1.0);

            const Point2 pmedian(pmedx, pmedy);

            //rotate(anticlockwise, el_median, el_end, getPosition);

            const auto[lpoly, rpoly] = bisect_polygon(domain, avg_vel, pmedian);

            bisected_parts.emplace_back(lpoly, el_begin,  el_median);
            bisected_parts.emplace_back(rpoly, el_median, el_end);
        }

        npart *= 2;
        partitions = bisected_parts;
    }

    return partitions;
}

}

}
#endif //NORCB_NORCB_HPP
