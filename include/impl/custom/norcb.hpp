//
// Created by xetql on 5/25/21.
//

#ifndef YALBB_NORCB_HPP
#define YALBB_NORCB_HPP
#include "impl/custom/geometry.hpp"
#include "geometry_utils.hpp"

#include <parallel/algorithm.hpp>
#include <tuple>
#include <iostream>
#include <mpi.h>
#include <vector>
#include <mpi.h>
#include <iterator>
#include <algorithm>
#include <utility>
namespace norcb {
template<class Real, class ForwardIt, class GetVelocityFunc>
Vector parallel_compute_average_velocity(ForwardIt el_begin, ForwardIt el_end, unsigned longest_axis, MPI_Comm comm,
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

struct NoRCB {

    int world_size, rank;
    Polygon domain;
    std::vector<Polygon> subdomains{};
    MPI_Comm comm;

    NoRCB(const Polygon& domain, MPI_Comm comm) : domain(domain), comm(comm) {
        MPI_Comm_size(comm, &world_size);
        MPI_Comm_rank(comm, &rank);
        subdomains.resize(world_size);
    }

    NoRCB(const Polygon& domain, std::vector<Polygon> subdomains, MPI_Comm comm) : domain(domain), comm(comm), subdomains(std::move(subdomains)) {
        MPI_Comm_size(comm, &world_size);
        MPI_Comm_rank(comm, &rank);
    }

    template<class Real>
    void get_owner(Real x, Real y, int* owner) {
        Point pt{x, y};

        const auto c = pt;
        for((*owner) = 0; (*owner) < world_size; (*owner)++){
            const auto& subdomain = this->subdomains.at(*owner);
            int belongs = 1;
            for(auto beg = subdomain.edges_begin(); beg != subdomain.edges_end(); beg++) {
                const auto& s = *beg;

                const auto&[a, b] = s;
                const auto ab = b-a;
                const auto ac = c-a;
                auto cross = (ab.x*ac.y - ab.y*ac.x);
                const auto sign  = std::signbit(cross);
                const auto is_on_boundary = std::fabs(cross) < acceptable_error;
                const auto side  = is_on_boundary || !sign;

                belongs &= side;
            }

            if(belongs) return;
        }
        *owner = -1;
    }

    template<class Real>
    void get_neighbors(Real x, Real y, Real z, double rc, int* PEs, int* num_found) const {
        Point p{x, y};
        const double sqrc = rc*rc;
        *num_found = 0;
        int i = 0;
        for(const Polygon& poly : subdomains) {
            for(auto beg = poly.edges_begin(); beg != poly.edges_end(); beg++) {
                if(sqr_dist(*beg, p) < sqrc) {
                    PEs[*num_found] = i;
                    *num_found = *num_found + 1;
                    break;
                }
            }
            i++;
        }
    }

};

template<class Real, class ForwardIt, class GetPositionFunc, class GetVelocityFunc>
void partition(NoRCB* lb_struct, unsigned P, ForwardIt el_begin, ForwardIt el_end,
               MPI_Datatype datatype, MPI_Comm comm,
          GetPositionFunc getPosition, GetVelocityFunc getVelocity) {

    std::vector<std::tuple<Polygon, ForwardIt, ForwardIt>> partitions {};

    partitions.emplace_back(lb_struct->domain, el_begin, el_end);

    unsigned npart = partitions.size();

    while (npart != P) {
        decltype(partitions) bisected_parts{};
        for (auto& partition : partitions) {
            auto&[domain, el_begin, el_end] = partition;
            const unsigned elements_in_subdomain = std::distance(el_begin, el_end);
            std::array<Real, 2> min{std::numeric_limits<Real>::max(), std::numeric_limits<Real>::max()};
            std::array<Real, 2> max{std::numeric_limits<Real>::lowest(), std::numeric_limits<Real>::lowest()};

            std::for_each(el_begin, el_end, [&min, &max, &getPosition] (const auto& e1) {
                auto position = getPosition(&e1);
                min.at(0) = std::min(position->at(0), min.at(0));
                max.at(0) = std::max(position->at(0), max.at(0));
                min.at(1) = std::min(position->at(1), min.at(1));
                max.at(1) = std::max(position->at(1), max.at(1));
            });

            MPI_Allreduce(MPI_IN_PLACE, min.data(), min.size(), par::get_mpi_type<typename decltype(min)::value_type>(), MPI_MIN, comm);
            MPI_Allreduce(MPI_IN_PLACE, max.data(), max.size(), par::get_mpi_type<typename decltype(max)::value_type>(),  MPI_MAX, comm);

            const unsigned target_axis = (max.at(0) - min.at(0)) < (max.at(1) - min.at(1));

            const Real ox = 0.0;
            const Real oy = 1.0;

            auto [avg_velx, avg_vely] = par_get_average_velocity<Real>(el_begin, el_end, target_axis, comm, getVelocity);

            auto theta   = get_angle(avg_velx, avg_vely, ox, oy);

            const auto clockwise     = get_rotation_matrix(theta);

            const auto anticlockwise = get_rotation_matrix(-theta);

            rotate(clockwise, el_begin, el_end, getPosition);

            const auto norm = std::sqrt(avg_velx*avg_velx+avg_vely*avg_vely);
            avg_velx = avg_velx / norm;
            avg_vely = avg_vely / norm;

            Real median;
            {
                std::remove_const_t<decltype(elements_in_subdomain)>  size {0};
                MPI_Allreduce(&elements_in_subdomain, &size, 1, par::get_mpi_type<decltype(size)>(), MPI_SUM, comm);
                auto midIdx = size / 2 - 1;
                auto el_median = par::find_nth(el_begin, el_end, midIdx, datatype, comm, [&getPosition](const auto& a, const auto& b){
                    auto posa = getPosition(&a);
                    auto posb = getPosition(&b);
                    return posa->at(0) < posb->at(0);
                });
                median = getPosition(&el_median)->at(0);
            }

            const auto c = std::count_if(el_begin, el_end, [&getPosition, &median](const auto& v) {
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

            const auto[lpoly, rpoly] = bisect_polygon(domain, avg_velx, avg_vely, pmedx, pmedy);

            bisected_parts.emplace_back(lpoly, el_begin,  el_median);
            bisected_parts.emplace_back(rpoly, el_median, el_end);
        }

        npart *= 2;
        partitions = bisected_parts;
    }

    for(auto i = 0; i < lb_struct->world_size; ++i) {
        auto& p = partitions.at(i);
        lb_struct->subdomains.at(i) = std::get<0>(p);
    }
}
}

#endif //YALBB_NORCB_HPP
