//
// Created by xetql on 3/17/21.
//

#ifndef NORCB_UTILS_HPP
#define NORCB_UTILS_HPP

template<class T>
constexpr T pi(){
    return T(3.14);
}

template<class T>
typename std::enable_if<!std::numeric_limits<T>::is_integer, bool>::type almost_equal(T x, T y, int ulp) {
    // the machine epsilon has to be scaled to the magnitude of the values used
    // and multiplied by the desired precision in ULPs (units in the last place)
    return std::fabs(x - y) <= std::numeric_limits<T>::epsilon() * std::fabs(x + y) * ulp
           // unless the result is subnormal
           || std::fabs(x - y) < std::numeric_limits<T>::min();
}

template<class Real>
int sign(Real x, int ulp = 2) {
    //return x < static_cast<Real>(0) ? static_cast<Real>(-1) : (almost_equal(x, static_cast<Real>(0), 2) ? static_cast<Real>(0) : static_cast<Real>(1));
    return (almost_equal(x, static_cast<Real>(0), ulp) ? static_cast<Real>(0) :
            x < static_cast<Real>(0) ? static_cast<Real>(-1) : static_cast<Real>(1));

}

#endif //NORCB_UTILS_HPP
