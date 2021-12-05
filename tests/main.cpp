//
// Created by xetql on 12/27/20.
//
#include <iostream>
#include "norcb.hpp"
#include <random>
#include "norcb.hpp"
#include <parallel/algorithm.hpp>
#define CGAL_DISABLE_ROUNDING_MATH_CHECK ON
using Real = float;
using Index= size_t;
using namespace norcb;
    template<unsigned N>
    struct Element {
        static const auto dimension = N;

        Index gid;
        Index lid;
        std::array<Real, N> position, velocity;

        constexpr Element(std::array<Real, N> p, std::array<Real, N> v, const Index gid, const Index lid) : gid(gid),
                                                                                                            lid(lid),
                                                                                                            position(p),
                                                                                                            velocity(v) {}

        constexpr Element() : gid(0), lid(0), position(), velocity() {}

        friend std::ostream &operator<<(std::ostream &os, const Element &element) {
            os << element.position.at(0);
            for (int i = 1; i < N; i++) {
                os << "," << element.position.at(i);
            }
            os << "; ";
            os << element.position.at(0);
            for (int i = 1; i < N; i++) {
                os << "," << element.velocity.at(i);
            }
            os << element.gid << ";" << element.lid;
            return os;
        }

        inline static MPI_Datatype register_datatype() {
            constexpr const bool UseDoublePrecision = std::is_same<Real, double>::value;
            MPI_Datatype element_datatype, vec_datatype, oldtype_element[2];

            MPI_Aint offset[2], lb, intex;

            int blockcount_element[2];

            // register particle element type
            constexpr int array_size = N;
            auto mpi_raw_datatype = UseDoublePrecision ? MPI_DOUBLE : MPI_FLOAT;

            MPI_Type_contiguous(array_size, mpi_raw_datatype, &vec_datatype);

            MPI_Type_commit(&vec_datatype);

            blockcount_element[0] = 2; //gid, lid
            blockcount_element[1] = 2; //position, velocity

            oldtype_element[0] = MPI_LONG_LONG;
            oldtype_element[1] = vec_datatype;

            MPI_Type_get_extent(MPI_LONG_LONG, &lb, &intex);

            offset[0] = static_cast<MPI_Aint>(0);
            offset[1] = blockcount_element[0] * intex;

            MPI_Type_create_struct(2, blockcount_element, offset, oldtype_element, &element_datatype);

            MPI_Type_commit(&element_datatype);

            return element_datatype;
        }

        inline static std::array<Real, N> *getElementPositionPtr(Element<N> *e) {
            return &(e->position);
        }

        inline static std::array<Real, N> *getElementVelocityPtr(Element<N> *e) {
            return &(e->velocity);
        }
    };

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    auto datatype = Element<2>::register_datatype();
    int wsize, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &wsize);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    auto d = init_domain(0., 0., 10., 10.);

    unsigned N = 10 + rank;

    std::random_device rd{};
    std::mt19937 gen(time(NULL));
    std::uniform_real_distribution<Real> dist(0.0, 1.0);

    std::vector<Real> x, y;
    std::vector<Real> vx, vy;
    std::vector<Element<2>> particles {};

    if(!rank)
        for(unsigned i = 0; i < 10000*wsize*N; ++i){
            std::array<Real, 2> pos = {dist(gen), dist(gen)};
            std::array<Real, 2> vel = {dist(gen), dist(gen)};
            particles.emplace_back(pos, vel, 0, 0);
        }

    CGAL::set_pretty_mode(std::cout);

    do_partition<Element<2>, Real>();

    MPI_Finalize();

    return EXIT_SUCCESS;
}