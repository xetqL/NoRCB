//
// Created by xetql on 5/25/21.
//

#ifndef YALBB_NORCB_HPP
#define YALBB_NORCB_HPP

#include "impl/CGAL/geometry.hpp"

#include <parallel/algorithm.hpp>
#include "geometry_utils.hpp"
#include <tuple>
#include <iostream>
#include <mpi.h>
#include <vector>
#include <mpi.h>

struct NoRCB {
    K trait {};
    int world_size, rank;
    Polygon2 domain;
    std::vector<Polygon2> subdomains{};
    MPI_Comm comm;

    NoRCB(const Polygon2& domain, MPI_Comm comm) : domain(domain), comm(comm) {
        MPI_Comm_size(comm, &world_size);
        MPI_Comm_rank(comm, &rank);
        subdomains.resize(world_size);
    }

    NoRCB(const Polygon2& domain, std::vector<Polygon2> subdomains, MPI_Comm comm) : domain(domain), comm(comm), subdomains(std::move(subdomains)) {
        MPI_Comm_size(comm, &world_size);
        MPI_Comm_rank(comm, &rank);
    }

    template<class Real>
    void get_owner(Real x, Real y, int* owner) {
        Point2 c(x,y);
        for((*owner) = 0; (*owner) < world_size; (*owner)++){
            const auto& subdomain = this->subdomains.at(*owner);
            int belongs= 1;
            for(auto beg = subdomain.edges_begin(); beg != subdomain.edges_end(); beg++) {
                Segment2 s(*beg);
                const auto& a = s.source();
                const auto& b = s.target();
                const auto ab = b-a;
                const auto ac = c-a;
                const auto cross = CGAL::to_double(ab.x()*ac.y() - ab.y()*ac.x());
                const auto sign  = std::signbit(cross);
                const auto is_on_boundary = std::fabs(cross) < acceptable_error;
                const auto side  = is_on_boundary || !sign;
                belongs &= side;
            }

            if(belongs) return;
        }

        throw std::logic_error(fmt("Point(%f, %f) belongs to no one", x, y));
    }

    template<class Real>
    void get_intersecting_domains(Real x1, Real x2, Real y1, Real y2, Real z1, Real z2, int* PEs, int* num_found) {
        Polygon2 interaction_square = init_domain(x1, y1, x2, y2);
        *num_found = 0;
        for(auto PE = 0; PE < world_size; ++PE){
            const Polygon2& subdomain = subdomains.at(PE);
            if(CGAL::do_intersect(subdomain, interaction_square)){
                PEs[*num_found] = PE;
                *num_found += 1;
            }
        }
    }

    template<class Real>
    void get_neighbors(Real x, Real y, Real z, double rc, int* PEs, int* num_found) const {
        Point2 p(x, y);
        const double sqrc = rc*rc;
        *num_found = 0;
        int i = 0;
        for(const Polygon2& poly : subdomains) {
            for(auto beg = poly.edges_begin(); beg != poly.edges_end(); beg++) {
                if(CGAL::squared_distance(*beg, p) < sqrc) {
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

    std::vector<std::tuple<Polygon2, ForwardIt, ForwardIt>> partitions {};

    partitions.emplace_back(lb_struct->domain, el_begin, el_end);

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
            //auto avg_vel = parallel_compute_average_velocity<Real>(el_begin, el_end, target_axis, comm, getPosition, getVelocity);
            auto [avg_x, avg_y] = par_get_average_velocity<Real>(el_begin, el_end, target_axis, comm, getPosition, getVelocity);
            Vector2 avg_vel(avg_x, avg_y);
            auto theta         = get_angle(avg_vel, origin);

            const auto clockwise     = get_rotation_matrix(theta);
            const auto anticlockwise = get_rotation_matrix(-theta);

            rotate(clockwise, el_begin, el_end, getPosition);

            const auto norm = std::sqrt(CGAL::to_double(avg_vel.x() * avg_vel.x() + avg_vel.y() * avg_vel.y()));
            avg_vel = avg_vel / norm;

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
    for(auto i = 0; i < lb_struct->world_size; ++i) {
        auto& p = partitions.at(i);
        lb_struct->subdomains.at(i) = std::get<0>(p);
    }

}
#endif //YALBB_NORCB_HPP
