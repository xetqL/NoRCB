//
// Created by xetql on 12/27/20.
//
#include <iostream>
#include "norcb.hpp"
#include <random>

#define CGAL_DISABLE_ROUNDING_MATH_CHECK ON
using namespace norcb;

struct Particle2 {
    using value_type = double;
    std::array<value_type, 2> position {}, velocity {};
    Particle2(std::array<value_type, 2> position, std::array<value_type, 2> velocity) : position(position),
        velocity(velocity) {}
};

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    MPI_Aint displacements[2]  = {offsetof(Particle2, position), offsetof(Particle2, velocity)};
    int block_lengths[2]  = {2, 2};
    MPI_Datatype types[2] = {MPI_DOUBLE, MPI_DOUBLE};
    MPI_Datatype datatype;
    MPI_Type_create_struct(2, block_lengths, displacements, types, &datatype);
    MPI_Type_commit(&datatype);
    int wsize, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &wsize);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    auto d = init_domain(0., 0., 1., 1.);

    unsigned N = 10 + rank;

    std::random_device rd{};
    std::mt19937 gen(rank);
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    std::vector<double> x, y;
    std::vector<double> vx, vy;
    std::vector<Particle2> particles {};

    for(unsigned i = 0; i < N; ++i){
        particles.push_back(Particle2({dist(gen),dist(gen)},{dist(gen),dist(gen)}));
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

    MPI_Finalize();

    return EXIT_SUCCESS;
}