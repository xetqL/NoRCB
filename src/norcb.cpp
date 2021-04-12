//
// Created by xetql on 12/27/20.
//

#include "norcb.hpp"

namespace norcb {

NoRCB* allocate_from(NoRCB* from) {
    return new NoRCB(from->subdomains, from->comm);
}
void destroy(NoRCB* lb) {
    delete lb;
}

namespace seq {
namespace{
Vector2 compute_average_velocity(std::vector<double> vx, std::vector<double> vy, bool axis) {
    return norcb::compute_average_velocity(std::move(vx), std::move(vy), axis, [](auto beg, auto end) {
        return std::make_tuple(std::accumulate(beg, end, 0.), std::distance(beg, end));
    });
}
}

std::vector<std::tuple<Polygon2,
        std::vector<double>, std::vector<double>,
        std::vector<double>, std::vector<double>>> partition(Integer P, std::vector<double> &x, std::vector<double> &y,
                                                             std::vector<double> &vx, std::vector<double> &vy,
                                                             const Polygon2 &domain) {
    std::vector<std::tuple<Polygon2,
            std::vector<double>, std::vector<double>,
            std::vector<double>, std::vector<double>>> partitions{};

    partitions.emplace_back(domain, std::move(x), std::move(y), std::move(vx), std::move(vy));
    while (P > 1) {
        decltype(partitions) bisected_parts{};
        for (auto &partition : partitions) {
            auto&[domain, x, y, vx, vy] = partition;
            const auto part_in_subdomain = x.size();

            const double minx = *std::min_element(x.cbegin(), x.cend());
            const double maxx = *std::max_element(x.cbegin(), x.cend());
            const double miny = *std::min_element(y.cbegin(), y.cend());
            const double maxy = *std::max_element(y.cbegin(), y.cend());

            const unsigned target_axis = ((maxx - minx) > (maxy - miny)) ? 0 : 1;

            const Vector2 origin(0., 1.);

            auto avg_vel = compute_average_velocity(vx, vy, target_axis);

            const auto theta = get_angle(avg_vel, origin);
            const auto clockwise = get_rotation_matrix(theta);
            const auto anticlockwise = get_rotation_matrix(-theta);

            rotate(clockwise, x, y);
            avg_vel = Vector2(
                    clockwise[0][0] * avg_vel.x() + clockwise[0][1] * avg_vel.y(),
                    clockwise[1][0] * avg_vel.x() + clockwise[1][1] * avg_vel.y()
            );

            const double norm = std::sqrt(avg_vel.x() * avg_vel.x() + avg_vel.y() * avg_vel.y());
            avg_vel = avg_vel / norm;

            double median;
            {
                std::vector<double> rx(x.begin(), x.end());
                unsigned midIdx = (rx.size() / 2);
                std::nth_element(rx.begin(), (rx.begin() + midIdx), rx.end());
                double a = *(rx.begin() + midIdx);
                if (rx.size() % 2) {
                    median = a;
                } else {
                    std::nth_element(rx.begin(), (rx.begin() + (midIdx - 1)), rx.end());
                    median = (a + *(rx.begin() + (midIdx - 1))) / 2;
                }
            }

            const auto c = std::count_if(x.begin(), x.end(), [median](auto v) { return v <= median; });

            std::vector<double> lx(c), ly(c), rx(part_in_subdomain - c), ry(part_in_subdomain - c);
            std::vector<double> lvx(c), lvy(c), rvx(part_in_subdomain - c), rvy(part_in_subdomain - c);

            unsigned l = 0, r = 0;
            for (unsigned i = 0; i < part_in_subdomain; ++i) {
                if (x[i] <= median) {
                    lx[l] = x[i];
                    ly[l] = y[i];
                    lvx[l] = vx[i];
                    lvy[l] = vy[i];
                    l++;
                } else {
                    rx[r] = x[i];
                    ry[r] = y[i];
                    rvx[r] = vx[i];
                    rvy[r] = vy[i];
                    r++;
                }
            }

            rotate(anticlockwise, lx, ly);

            const auto[pmedx, pmedy] = rotate(anticlockwise, median, 1.0);

            const Point2 pmedian(pmedx, pmedy);

            rotate(anticlockwise, rx, ry);

            const auto[lpoly, rpoly] = bisect_polygon(domain, avg_vel, pmedian);

            bisected_parts.emplace_back(lpoly, lx, ly, lvx, lvy);
            bisected_parts.emplace_back(rpoly, rx, ry, rvx, rvy);

        }
        P /= 2;
        partitions = bisected_parts;
    }

    return partitions;
}

}

namespace parallel {
namespace {
Vector2 compute_average_velocity(std::vector<double> vx, std::vector<double> vy, bool axis, MPI_Comm comm) {
    return norcb::compute_average_velocity(std::move(vx), std::move(vy), axis, [comm](auto beg, auto end) {
        auto acc = std::accumulate(beg, end, 0.);
        MPI_Allreduce(MPI_IN_PLACE, &acc, 1, par::get_mpi_type<decltype(acc)>(), MPI_SUM, comm);
        auto size = std::distance(beg, end);
        MPI_Allreduce(MPI_IN_PLACE, &size, 1, par::get_mpi_type<decltype(size)>(), MPI_SUM, comm);
        return std::make_tuple(acc, size);
    });
}
}
std::vector<std::tuple<Polygon2,
    std::vector<double>, std::vector<double>,
    std::vector<double>, std::vector<double>>> partition(Integer P,
     std::vector<double> &x,  std::vector<double> &y,
     std::vector<double> &vx, std::vector<double> &vy,
     const Polygon2 &domain, MPI_Comm comm) {
    std::vector<std::tuple<Polygon2,
            std::vector<double>, std::vector<double>,
            std::vector<double>, std::vector<double>>> partitions {};

    partitions.emplace_back(domain, std::move(x), std::move(y), std::move(vx), std::move(vy));
    unsigned npart = partitions.size();
    while (npart != P) {
        decltype(partitions) bisected_parts{};
        for (auto &partition : partitions) {
            auto&[domain, x, y, vx, vy] = partition;
            const auto part_in_subdomain = x.size();

            const double minx = *std::min_element(x.cbegin(), x.cend());
            const double maxx = *std::max_element(x.cbegin(), x.cend());
            const double miny = *std::min_element(y.cbegin(), y.cend());
            const double maxy = *std::max_element(y.cbegin(), y.cend());

            const unsigned target_axis = ((maxx - minx) > (maxy - miny)) ? 0 : 1;

            const Vector2 origin(0., 1.);

            auto avg_vel = compute_average_velocity(vx, vy, target_axis, comm);

            const auto theta         = get_angle(avg_vel, origin);
            const auto clockwise     = get_rotation_matrix(theta);
            const auto anticlockwise = get_rotation_matrix(-theta);

            rotate(clockwise, x, y);

            const double norm = std::sqrt(avg_vel.x() * avg_vel.x() + avg_vel.y() * avg_vel.y());
            avg_vel = avg_vel / norm;

            double median;

            {
                std::vector<double> rx(x.begin(), x.end());
                auto size = rx.size();
                MPI_Allreduce(MPI_IN_PLACE, &size, 1, par::get_mpi_type<decltype(size)>(), MPI_SUM, comm);
                unsigned midIdx = size / 2;
                median = par::find_nth(rx.begin(), rx.end(), midIdx, comm);
            }

            const auto c = std::count_if(x.begin(), x.end(), [median](auto v) { return v <= median; });

            std::vector<double> lx(c), ly(c), rx(part_in_subdomain - c), ry(part_in_subdomain - c);
            std::vector<double> lvx(c), lvy(c), rvx(part_in_subdomain - c), rvy(part_in_subdomain - c);

            unsigned l = 0, r = 0;
            for (unsigned i = 0; i < part_in_subdomain; ++i) {
                if (x[i] <= median) {
                    lx[l] = x[i];
                    ly[l] = y[i];
                    lvx[l] = vx[i];
                    lvy[l] = vy[i];
                    l++;
                } else {
                    rx[r] = x[i];
                    ry[r] = y[i];
                    rvx[r] = vx[i];
                    rvy[r] = vy[i];
                    r++;
                }
            }

            rotate(anticlockwise, lx, ly);

            const auto[pmedx, pmedy] = rotate(anticlockwise, median, 1.0);

            const Point2 pmedian(pmedx, pmedy);
            rotate(anticlockwise, rx, ry);

            const auto[lpoly, rpoly] = bisect_polygon(domain, avg_vel, pmedian);

            bisected_parts.emplace_back(lpoly, lx, ly, lvx, lvy);
            bisected_parts.emplace_back(rpoly, rx, ry, rvx, rvy);

        }
        npart *= 2;
        partitions = bisected_parts;
    }

    return partitions;
}

}
}