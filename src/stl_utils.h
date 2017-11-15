#ifndef STL_UTILS_H
#define STL_UTILS_H

#include <array>
#include <bitset>
#include <iostream>
#include <vector>

template <class T, class F>
void rearranging_erase_if(std::vector<T> &c, F f)
{
    for (size_t i = c.size(); i--;) {
        if (f(c[i])) {
            c[i] = c.back();
            c.pop_back();
        }
    }
}

template <size_t N>
bool operator<(const std::bitset<N> &lhs, const std::bitset<N> &rhs)
{
    return (lhs | rhs) == rhs;
}

template <class T>
std::ostream &operator<<(std::ostream &os, const std::vector<T> &rhs)
{
    os << "{";
    if (rhs.size() != 0)
        os << rhs[0];
    for (size_t i = 1; i < rhs.size(); ++i)
        os << ", " << rhs[i];
    return os << "}";
}

template <class T1, class T2>
std::ostream &operator<<(std::ostream &os, const std::pair<T1, T2> &rhs)
{
    return os << "(" << rhs.first << ", " << rhs.second << ")";
}

#endif
