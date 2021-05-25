//
// Created by xetql on 12/27/20.
//

#ifndef NORCB_NORCB_HPP
#define NORCB_NORCB_HPP

#ifdef USE_CGAL
#include "impl/CGAL/norcb.hpp"
#else
#include "impl/custom/norcb.hpp"
#endif

namespace norcb {

inline NoRCB* allocate_from(NoRCB* from) {
    return new NoRCB(from->domain, from->subdomains, from->comm);
}
inline void destroy(NoRCB* lb) {
    delete lb;
}

}
#endif //NORCB_NORCB_HPP
