//
// Created by xetql on 12/27/20.
//
#include <iostream>
#include "norcb.hpp"
int main(int argc, char** argv) {
    auto d = init_domain(0., 0., 1., 1.);
    auto [p1,p2] = bisect_polygon(d, Point2(1., 1.), Vector2(1., 1.));
    std::cout << p1 << std::endl;
    std::cout << p2 << std::endl;
    return EXIT_SUCCESS;
}