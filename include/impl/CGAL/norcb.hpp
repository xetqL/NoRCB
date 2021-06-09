//
// Created by xetql on 5/25/21.
//

#ifndef YALBB_NORCB_HPP
#define YALBB_NORCB_HPP

#include "impl/CGAL/geometry.hpp"

#include <parallel/algorithm.hpp>
#include <parallel/blocklb.hpp>
#include <parallel/geomedian.hpp>
#include "geometry_utils.hpp"
#include <tuple>
#include <iostream>
#include <mpi.h>
#include <vector>
#include <mpi.h>
#include <any>

namespace norcb {

struct BisectionTree {
    std::any median;
    BisectionTree *left = nullptr,
            *right = nullptr;

    explicit BisectionTree(std::any val): median(val){}

    template<class T>
    void addLeft(T val) {
        left = new BisectionTree(std::move(val));
    }
    template<class T>
    void addRight(T val){
        right = new BisectionTree(std::move(val));
    }

    ~BisectionTree() {
        delete left;
        delete right;
    }
};

struct NoRCB {
    int world_size, rank;
    Polygon2 domain;
    std::vector<Polygon2> subdomains{};
    MPI_Comm comm;

    BisectionTree* bisections = nullptr;

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
    using T = typename ForwardIt::value_type;
    int rank, nprocs;
    MPI_Comm_size(comm, &nprocs);
    MPI_Comm_rank(comm, &rank);
    std::vector<std::tuple<Polygon2, ForwardIt, ForwardIt>> partitions {};
    std::vector<BisectionTree**> bisection_queue = {&lb_struct->bisections};
    partitions.emplace_back(lb_struct->domain, el_begin, el_end);

    unsigned npart = 1;
    while (npart != P) {
        decltype(partitions) bisected_parts{};
        decltype(bisection_queue) next_bisection_ptr{};
        for (auto ipart = 0; ipart < npart; ++ipart) {
            auto &partition = partitions.at(ipart);
            auto ptr_bisection = bisection_queue.at(ipart);

            auto&[domain, _el_begin, _el_end] = partition;

            std::vector<T> current_data(_el_begin, _el_end);
            blocklb(rank, nprocs, current_data, datatype, comm);
            auto el_begin = current_data.begin();
            auto el_end   = current_data.end();

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
            auto [avg_x, avg_y] = par_get_average_velocity<Real>(el_begin, el_end, target_axis, comm, getVelocity);
            Vector2 avg_vel(avg_x, avg_y);
            auto theta               = get_angle(avg_vel, origin);
            const auto clockwise     = get_rotation_matrix(theta);
            const auto anticlockwise = get_rotation_matrix(-theta);

            rotate(clockwise, el_begin, el_end, getPosition);

            const auto norm = std::sqrt(CGAL::to_double(avg_vel.x() * avg_vel.x() + avg_vel.y() * avg_vel.y()));
            avg_vel = avg_vel / norm;

            Real median;
            std::optional<Real> pivot_hint;

            if(*ptr_bisection) {
                pivot_hint = std::any_cast<Real>((*ptr_bisection)->median);
            } else {
                pivot_hint = std::nullopt;
            }

            {
                std::remove_const_t<decltype(elements_in_subdomain)>  size {0};
                MPI_Allreduce(&elements_in_subdomain, &size, 1, par::get_mpi_type<decltype(size)>(), MPI_SUM, comm);
                auto midIdx = size / 2 - 1;

                median = par::find_spatial_median(rank, nprocs, el_begin, el_end, 0.001, comm,
                                                  [getPosition](const auto& v){return getPosition(&v)->at(0);},
                                                  std::nullopt);

                if(!(*ptr_bisection)) {
                    *ptr_bisection = new BisectionTree(median);
                } else {
                    (*ptr_bisection)->median = median;
                }
            }

            rotate(clockwise, _el_begin, _el_end, getPosition);

            auto el_median_it = std::partition(_el_begin, _el_end, [&getPosition, &median](const auto& el){
                return getPosition(&el)->at(0) <= median;
            });

            rotate(anticlockwise, _el_begin, _el_end, getPosition);

            const auto[pmedx, pmedy] = rotate(anticlockwise, median, 1.0);

            const Point2 pmedian(pmedx, pmedy);

            const auto[lpoly, rpoly] = bisect_polygon(domain, avg_vel, pmedian);

            bisected_parts.emplace_back(lpoly, _el_begin,  el_median_it);
            next_bisection_ptr.push_back(&(*ptr_bisection)->left);
            bisected_parts.emplace_back(rpoly, el_median_it, _el_end);
            next_bisection_ptr.push_back(&(*ptr_bisection)->right);
        }

        npart *= 2;
        partitions = bisected_parts;
        bisection_queue = next_bisection_ptr;
    }
    for(auto i = 0; i < lb_struct->world_size; ++i) {
        auto& p = partitions.at(i);
        lb_struct->subdomains.at(i) = std::get<0>(p);
    }

}
}
#endif //YALBB_NORCB_HPP
