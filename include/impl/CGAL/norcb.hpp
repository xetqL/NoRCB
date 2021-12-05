//
// Created by xetql on 5/25/21.
//

#ifndef YALBB_NORCB_HPP
#define YALBB_NORCB_HPP

#include "impl/CGAL/geometry.hpp"

#include <parallel/algorithm.hpp>
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

template<class Real, class RandomIt, class GetPositionFunc, class GetVelocityFunc>
void partition(NoRCB* lb_struct, unsigned P, RandomIt el_begin, RandomIt el_end,
               MPI_Datatype datatype, MPI_Comm comm,
               GetPositionFunc getPosition, GetVelocityFunc getVelocity) {
    using T = typename RandomIt::value_type;
    int rank, nprocs;
    MPI_Comm_size(comm, &nprocs);
    MPI_Comm_rank(comm, &rank);
    std::vector<std::tuple<Polygon2, RandomIt, RandomIt>> partitions {};
    std::vector<BisectionTree**> bisection_queue = {&lb_struct->bisections};
    partitions.emplace_back(lb_struct->domain, el_begin, el_end);

    unsigned npart = 1;
    int iter = 0;
    while (npart != P) {
        decltype(partitions) bisected_parts{};
        decltype(bisection_queue) next_bisection_ptr{};
        for (auto ipart = 0; ipart < npart; ++ipart) {
            auto &partition = partitions.at(ipart);
            auto ptr_bisection = bisection_queue.at(ipart);

            auto&[domain, el_begin, el_end] = partition;

            unsigned elements_in_subdomain = std::distance(el_begin, el_end);

            // compute minx, maxx, miny, maxy

            std::array<Real, 2> mins = {std::numeric_limits<Real>::max(), std::numeric_limits<Real>::max()};
            std::array<Real, 2> maxs = {std::numeric_limits<Real>::lowest(), std::numeric_limits<Real>::lowest()};

            std::for_each(el_begin, el_end, [&mins, &maxs, getPosition] (const auto& e1) {
                auto position = getPosition(&e1);
                mins[0] = std::min(position->at(0), mins[0]);
                maxs[0] = std::max(position->at(0), maxs[0]);
                mins[1] = std::min(position->at(1), mins[1]);
                maxs[1] = std::max(position->at(1), maxs[1]);
            });

            MPI_Allreduce(MPI_IN_PLACE, mins.data(), 2, par::get_mpi_type<Real>(), MPI_MIN, comm);
            MPI_Allreduce(MPI_IN_PLACE, maxs.data(), 2, par::get_mpi_type<Real>(), MPI_MAX, comm);

            // compute longest axis
            unsigned longest_axis = ((maxs[0] - mins[0]) < (maxs[1] - mins[1])) ? 1 : 0;

            const Vector2 origin(0.0, 1.0);

            // compute average velocity vector
            auto [avg_x, avg_y] = par_get_average_velocity<Real>(el_begin, el_end, longest_axis, comm, getVelocity);

            auto norm = std::sqrt( CGAL::to_double(avg_x*avg_x + avg_y*avg_y) );
            Vector2 avg_vel(avg_x, avg_y);

            if (norm < 1e-3) {
                avg_vel = Vector2(longest_axis, (longest_axis + 1) % 2) ;
            }

            // get angle between origin and velocity vector
            auto theta               = get_angle(avg_vel, origin);
            // get rotation matrix clockwise and anticlockwise
            const auto clockwise     = get_rotation_matrix(theta);
            const auto anticlockwise = get_rotation_matrix(-theta);

            // rotate data such that they are perpendicular to the cut axis (Y-axis)
            rotate(clockwise, el_begin, el_end, getPosition);

            // normalize velocity vector
            norm = std::sqrt(CGAL::to_double(avg_vel.x() * avg_vel.x() + avg_vel.y() * avg_vel.y()));
            avg_vel = avg_vel / norm;

            // compute median
            Real median;
            std::optional<Real> pivot_hint;
            if(*ptr_bisection) {
                pivot_hint = std::any_cast<Real>((*ptr_bisection)->median);
            } else {
                pivot_hint = std::nullopt;
            }

            RandomIt el_median_it;
            {
                decltype(elements_in_subdomain)  size;
                MPI_Allreduce(&elements_in_subdomain, &size, 1, par::get_mpi_type<decltype(size)>(), MPI_SUM, comm);

                auto midIdx = size / 2 - 1;
                auto el_median = par::find_nth(el_begin, el_end, midIdx, datatype, comm, [&getPosition](const auto& a, const auto& b){
                    auto posa = getPosition(&a);
                    auto posb = getPosition(&b);
                    return posa->at(0) < posb->at(0);
                });

                median = getPosition(&el_median)->at(0);
                // store median value for later use
                if(!(*ptr_bisection)) {
                    *ptr_bisection = new BisectionTree(median);
                } else {
                    (*ptr_bisection)->median = median;
                }
            } // end of median computation
            const auto c = std::count_if(el_begin, el_end, [&getPosition, &median](const auto& v) {
                    auto pos = getPosition(&v);
                    return pos->at(0) <= median;
                });
            el_median_it = el_begin + c;

            // rotate elements backwards
            rotate(anticlockwise, el_begin, el_end, getPosition);
            // rotate the median point backward
            const auto[pmedx, pmedy] = rotate(anticlockwise, median, 1.0);
            // Create the median point
            const Point2 pmedian(pmedx, pmedy);
            // bisect the current polygon using the median point and the velocity vector
            const auto[lpoly, rpoly] = bisect_polygon(domain, avg_vel, pmedian);
            // store the sub-domains for recursive partitioning
            bisected_parts.emplace_back(lpoly, el_begin,  el_median_it);
            next_bisection_ptr.push_back(&(*ptr_bisection)->left);
            bisected_parts.emplace_back(rpoly, el_median_it, el_end);
            next_bisection_ptr.push_back(&(*ptr_bisection)->right);
        }
        // number of partition increased by two
        npart *= 2;
        // sub-domains are the domain we will cut in the next iteration
        partitions = bisected_parts;
        bisection_queue = next_bisection_ptr;
        iter++;
    }

    // Update the load balancing structure
    for(auto i = 0; i < lb_struct->world_size; ++i) {
        auto& p = partitions.at(i);
        lb_struct->subdomains.at(i) = std::get<0>(p);
    }

}

inline bool is_set(int num, unsigned bit) {
    return (num>>bit) & 1;
}
inline int flip(int num, unsigned bit) {
    return num ^ (1 << bit);
}
template<class FwIt>
void sendrecv(std::vector<typename FwIt::value_type>& recv_buf, int* rcount, FwIt begin, FwIt end, MPI_Datatype datatype, int dest, MPI_Comm comm){
    const auto send_size = std::distance(begin, end);
    MPI_Request sreq{};
    MPI_Status status{};
    MPI_Isend(&(*begin), send_size, datatype, dest, 1223, comm, &sreq);
    MPI_Probe(dest, 1223, comm, &status);
    MPI_Get_count(&status, datatype, rcount);
    recv_buf.reserve(*rcount);
    MPI_Recv(recv_buf.data(), *rcount, datatype, dest, 1223, comm, MPI_STATUS_IGNORE);
    MPI_Wait(&sreq, MPI_STATUS_IGNORE);
}

template<class T, class Real, class GetPositionFunc, class GetVelocityFunc>
void do_partition(NoRCB* lb_struct, unsigned P, std::vector<T>& elements,
               MPI_Datatype datatype, MPI_Comm world_comm,
               GetPositionFunc getPosition, GetVelocityFunc getVelocity) {
    auto el_begin = elements.begin();
    auto el_end   = elements.end();
    unsigned elements_in_subdomain = std::distance(el_begin, el_end), total;
    MPI_Allreduce(&elements_in_subdomain, &total, 1, par::get_mpi_type<decltype(elements_in_subdomain)>(), MPI_SUM,  world_comm);

    std::vector<T> recv_buf;
    recv_buf.reserve(total / lb_struct->world_size);

    int rank, nprocs;

    MPI_Comm current_comm;
    MPI_Comm_dup(world_comm, &current_comm);

    auto domain = lb_struct->domain;

    auto iteration_index = 0u;

    auto n_partitions    = 1u;

    while (n_partitions < lb_struct->world_size)  {
        el_begin = elements.begin();
        el_end   = elements.end();
        MPI_Comm_rank(current_comm, &rank);
        MPI_Comm_size(current_comm, &nprocs);
        elements_in_subdomain = std::distance(el_begin, el_end);

        // compute minx, maxx, miny, maxy
        std::array<Real, 4> boundaries = {std::numeric_limits<Real>::max(),    std::numeric_limits<Real>::max(),
                                          std::numeric_limits<Real>::lowest(), std::numeric_limits<Real>::lowest()};

        std::for_each(el_begin, el_end, [&boundaries, getPosition] (const auto& e1) {
            auto position = getPosition(&e1);
            boundaries[0] = std::min(position->at(0), boundaries[0]);
            boundaries[1] = std::min(position->at(1), boundaries[1]);
            boundaries[2] = std::max(position->at(0), boundaries[2]);
            boundaries[3] = std::max(position->at(1), boundaries[3]);
        });

        std::array<Real, 2> mins = {std::numeric_limits<Real>::max(), std::numeric_limits<Real>::max()};
        std::array<Real, 2> maxs = {std::numeric_limits<Real>::lowest(), std::numeric_limits<Real>::lowest()};

        std::for_each(el_begin, el_end, [&mins, &maxs, getPosition] (const auto& e1) {
            auto position = getPosition(&e1);
            mins[0] = std::min(position->at(0), mins[0]);
            maxs[0] = std::max(position->at(0), maxs[0]);
            mins[1] = std::min(position->at(1), mins[1]);
            maxs[1] = std::max(position->at(1), maxs[1]);
        });

        MPI_Allreduce(MPI_IN_PLACE, mins.data(), 2, par::get_mpi_type<Real>(), MPI_MIN, world_comm);
        MPI_Allreduce(MPI_IN_PLACE, maxs.data(), 2, par::get_mpi_type<Real>(), MPI_MAX, world_comm);

        // compute longest axis
        unsigned longest_axis = ((maxs[0] - mins[0]) < (maxs[1] - mins[1])) ? 1 : 0;

        const Vector2 origin(0.0, 1.0);

        // compute average velocity vector
        auto [avg_x, avg_y] = par_get_average_velocity<Real>(el_begin, el_end, longest_axis, current_comm, getVelocity);

        auto norm = std::sqrt( CGAL::to_double(avg_x*avg_x + avg_y*avg_y) );
        Vector2 avg_vel(avg_x, avg_y);

        if (norm < 1e-3) {
            avg_vel = Vector2(longest_axis, (longest_axis + 1) % 2) ;
        }

        // get angle between origin and velocity vector
        auto theta                      = get_angle(avg_vel, origin);
        // get rotation matrix clockwise and anticlockwise
        const auto clockwise     = get_rotation_matrix(theta);
        const auto anticlockwise = get_rotation_matrix(-theta);

        // rotate data such that they are perpendicular to the cut axis (Y-axis)
        rotate(clockwise, el_begin, el_end, getPosition);

        // normalize velocity vector
        norm = std::sqrt(CGAL::to_double(avg_vel.x() * avg_vel.x() + avg_vel.y() * avg_vel.y()));
        avg_vel = avg_vel / norm;

        // compute median
        Real median;

        decltype(elements_in_subdomain)  size;
        MPI_Allreduce(&elements_in_subdomain, &size, 1, par::get_mpi_type<decltype(size)>(), MPI_SUM, current_comm);

        auto midIdx = size / 2 - 1;
        auto el_median = par::find_nth(el_begin, el_end, midIdx, datatype, current_comm, [&getPosition](const auto& a, const auto& b){
            auto posa = getPosition(&a);
            auto posb = getPosition(&b);
            return posa->at(0) < posb->at(0);
        });

        median = getPosition(&el_median)->at(0);

        const auto n_below_median = std::count_if(el_begin, el_end, [&getPosition, &median](const auto& v) {
            //auto pos = getPosition(&v);
            return getPosition(&v)->at(0) <= median;
        });

        auto el_median_it = el_begin + n_below_median;

        // rotate elements backwards
        rotate(anticlockwise, el_begin, el_end, getPosition);
        // rotate the median point backward
        const auto[pmedx, pmedy] = rotate(anticlockwise, median, 1.0);
        // Create the median point
        const Point2 pmedian(pmedx, pmedy);

        // bisect the current polygon using the median point and the velocity vector
        const auto[lpoly, rpoly] = bisect_polygon(domain, avg_vel, pmedian);

        decltype(el_begin) beg, end;
        const auto group = is_set(rank, 0);

        if (group) {
            domain = lpoly;
            // send these
            beg    = el_median_it;
            end    = el_end;
        } else {
            domain = rpoly;
            // send these
            beg    = el_begin;
            end    = el_median_it;
        }
        const int neighbor = flip(rank, 0);
        int count;
        sendrecv(recv_buf, &count, beg, end, datatype, neighbor, current_comm);

        elements.erase(beg, end);
        elements.reserve(elements.size() + count);

        std::move(recv_buf.data(), recv_buf.data()+count,  std::back_inserter(elements));

        iteration_index++;
        n_partitions *= 2;

        MPI_Comm split;
        MPI_Comm_split(current_comm, group, rank, &split);
        MPI_Comm_free(&current_comm);

        current_comm = split;

    }

    MPI_Comm_free(&current_comm);

    auto beg  = domain.vertices_begin();
    auto end  = domain.vertices_end();

    std::vector<Real> vertices(domain.size() * 2, 0.0);

    for(int i=0; i < domain.size(); beg++, i++) {
        const auto& p= *beg;
        vertices.at(2*i)   = p.x();
        vertices.at(2*i+1) = p.y();
    }

    std::vector<Real> rcv_vertices(lb_struct->world_size * vertices.size());
    MPI_Allgather(vertices.data(), vertices.size(), par::get_mpi_type<Real>(), rcv_vertices.data(), vertices.size(), par::get_mpi_type<Real>(), world_comm);
    for (int i = 0; i < lb_struct->world_size; ++i) {
        Polygon2 poly;
        for(int vpos = 0; vpos < vertices.size(); vpos+=2) {
            poly.push_back(
                    Point2(rcv_vertices.at(i * vertices.size() + vpos),
                           rcv_vertices.at(i * vertices.size() + vpos + 1)));
        }
        lb_struct->subdomains.at(i) = poly;
    }
}

}
#endif //YALBB_NORCB_HPP
