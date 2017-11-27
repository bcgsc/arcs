#ifndef _PAIR_HASH_H_
#define _PAIR_HASH_H_ 1

#include <functional>

/** Hash functor for std::pair of two arbitrary types */
struct PairHash {
    template <typename FirstT, typename SecondT>
    std::size_t operator()(const std::pair<FirstT, SecondT>& pair) const {
        return std::hash<FirstT>()(pair.first)
            ^ std::hash<SecondT>()(pair.second);
    }
};

#endif
