//
// Created by xetql on 12/27/20.
//
#include <iostream>
#include "norcb.hpp"
#include <random>
#include "norcb.hpp"
#define CGAL_DISABLE_ROUNDING_MATH_CHECK ON
using namespace norcb;
/*
struct Particle2 {
    using value_type = double;
    std::array<value_type, 2> position {}, velocity {};
    Particle2(std::array<value_type, 2> position, std::array<value_type, 2> velocity) : position(position),
        velocity(velocity) {}
};
*/
int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    auto datatype = elements::register_datatype<2>();
    int wsize, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &wsize);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    auto d = init_domain(0., 0., 10., 10.);

    unsigned N = 10 + rank;

    std::random_device rd{};
    std::mt19937 gen(time(NULL));
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    std::vector<double> x, y;
    std::vector<double> vx, vy;
    std::vector<elements::Element<2>> particles {};

    if(!rank)
        for(unsigned i = 0; i < 10000*wsize*N; ++i){
            std::array<double, 2> pos = {dist(gen), dist(gen)};
            std::array<double, 2> vel = {dist(gen), dist(gen)};
            particles.push_back(elements::Element<2>(pos, vel, 0, 0));
        }

    CGAL::set_pretty_mode(std::cout);

    auto parts = parallel::partition<double>(wsize, particles.begin(), particles.end(), d, datatype, MPI_COMM_WORLD, [](auto* p){
        return &p->position;
    }, [](auto* p){
        return &p->velocity;
    });

    std::vector<unsigned> domain_size(wsize, 0);
    std::for_each(parts.begin(), parts.end(),  [&domain_size, i=0] (const auto& p) mutable {
        domain_size.at(i) += (std::get<2>(p) - std::get<1>(p));
        i++;
    });

    MPI_Allreduce(MPI_IN_PLACE, domain_size.data(), wsize, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);

    if(!rank){
        std::copy(domain_size.begin(), domain_size.end(), std::ostream_iterator<int>(std::cout, " "));
        auto max = (double) *std::max_element(domain_size.begin(), domain_size.end());
        auto avg = (double) std::accumulate(domain_size.begin(), domain_size.end(), 0u) / domain_size.size();
        std::cout << ser::fmt("\nPercent Imbalance = %f\nImbalance Loss = %f", (max / avg)-1.0, max-avg) << std::endl;
    }


    particles.clear();
    if(rank % 2)
        for(unsigned i = 0; i < 5000*wsize*N; ++i) {
            std::array<double, 2> pos = {dist(gen), dist(gen)};
            std::array<double, 2> vel = {dist(gen), dist(gen)};
            particles.push_back(elements::Element<2>(pos, vel, 0, 0));
        }

    parts = parallel::partition<double>(wsize, particles.begin(), particles.end(), d, datatype, MPI_COMM_WORLD, [](auto* p){
        return &p->position;
    }, [](auto* p){
        return &p->velocity;
    });
    parts = parallel::partition<double>(wsize, particles.begin(), particles.end(), d, datatype, MPI_COMM_WORLD, [](auto* p){
        return &p->position;
    }, [](auto* p){
        return &p->velocity;
    });
    
    std::for_each(parts.begin(), parts.end(),  [&domain_size, i=0] (const auto& p) mutable {
        domain_size.at(i) += (std::get<2>(p) - std::get<1>(p));
        i++;
    });

    MPI_Allreduce(MPI_IN_PLACE, domain_size.data(), wsize, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);

    if(!rank){
        std::copy(domain_size.begin(), domain_size.end(), std::ostream_iterator<int>(std::cout, " "));
        auto max = (double) *std::max_element(domain_size.begin(), domain_size.end());
        auto avg = (double) std::accumulate(domain_size.begin(), domain_size.end(), 0u) / domain_size.size();
        std::cout << ser::fmt("\nPercent Imbalance = %f\nImbalance Loss = %f", (max / avg)-1.0, max-avg) << std::endl;
    }


    MPI_Finalize();

    return EXIT_SUCCESS;
}