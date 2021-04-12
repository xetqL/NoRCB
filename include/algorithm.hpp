//
// Created by xetql on 3/17/21.
//

#ifndef NORCB_ALGORITHM_HPP
#define NORCB_ALGORITHM_HPP

/**
 * Distinct without sorting worst is O(N^2) and best is O(N)
 * @tparam ForwardIt
 * @tparam BinaryComp
 * @param first
 * @param last
 * @param eq
 * @return
 */
template<class ForwardIt, class BinaryComp>
ForwardIt distinct(ForwardIt first, ForwardIt last, BinaryComp eq) {
    if (first == last) return last;
    while (first != last) {
        ForwardIt next = std::next(first);
        while (next != last) {
            if (eq(*first, *next)) { // not alone, swap and erase
                *next = *(last - 1);
                last--;
            } else next++;
        }
        first++;
    }
    return last;
}

#endif //NORCB_ALGORITHM_HPP
