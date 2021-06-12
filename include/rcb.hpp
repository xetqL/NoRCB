//
// Created by xetql on 6/10/21.
//

#pragma once
#ifdef USE_CGAL
#include "impl/CGAL/rcb.hpp"
#else
static_assert(false);
#endif
namespace rcb {

inline RCB* allocate_from(RCB* from) {
    return new RCB(from->domain, from->subdomains, from->comm);
}
inline void destroy(RCB* lb) {
    delete lb;
}

}