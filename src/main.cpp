//
// Created by xetql on 12/27/20.
//
#include <iostream>
#include "norcb.hpp"
#include <random>

#define CGAL_DISABLE_ROUNDING_MATH_CHECK ON
using namespace norcb;
int main(int argc, char** argv) {
    auto d = init_domain(0., 0., 1., 1.);

    std::random_device rd;

    std::uniform_real_distribution<double> dist(0.0, 1.0);

    std::vector<double> x,y;
    std::vector<double> vx,vy;

    for(unsigned i = 0; i < 8'010'0; ++i){
         x.push_back(dist(rd));
         y.push_back(dist(rd));
        vx.push_back(dist(rd));
        vy.push_back(dist(rd));
    }

    auto parts = seq::partition(8, x, y, vx, vy, d);

    std::for_each(parts.begin(), parts.end(), [](auto p){std::cout << std::get<1>(p).size() << " "; });

    return EXIT_SUCCESS;
}