#ifndef VECTOR_H
#define VECTOR_H

#include <array>
#include <cstddef>
#include <iostream>

#include "Global.h"

using Vector = std::array<Long, dim>;

inline Vector operator+(const Vector &lhs, const Vector &rhs)
{
    Vector ret{};
    for (size_t i = 0; i < dim; ++i)
        ret[i] = lhs[i] + rhs[i];
    return ret;
}

inline Vector operator-(const Vector &lhs, const Vector &rhs)
{
    Vector ret{};
    for (size_t i = 0; i < dim; ++i)
        ret[i] = lhs[i] - rhs[i];
    return ret;
}

inline Vector operator-(const Vector &rhs)
{
    Vector ret{};
    for (size_t i = 0; i < dim; ++i)
        ret[i] = -rhs[i];
    return ret;
}

inline Long operator*(const Vector &lhs, const Vector &rhs)
{
    Long ret{};
    for (size_t i = 0; i < dim; ++i)
        ret += lhs[i] * rhs[i];
    return ret;
}

inline Vector operator*(const Vector &lhs, Long rhs)
{
    Vector ret{};
    for (size_t i = 0; i < dim; ++i)
        ret[i] = lhs[i] * rhs;
    return ret;
}

inline Vector operator/(const Vector &lhs, Long rhs)
{
    Vector ret{};
    for (size_t i = 0; i < dim; ++i)
        ret[i] = lhs[i] / rhs;
    return ret;
}

inline Vector &operator/=(Vector &lhs, Long rhs)
{
    for (size_t i = 0; i < dim; ++i)
        lhs[i] /= rhs;
    return lhs;
}

inline Vector operator*(Long lhs, const Vector &rhs)
{
    Vector ret{};
    for (size_t i = 0; i < dim; ++i)
        ret[i] = rhs[i] * lhs;
    return ret;
}

inline std::ostream &operator<<(std::ostream &os, const Vector &rhs)
{
    os << "(";
    if (dim != 0)
        os << rhs[0];
    for (size_t i = 1; i < dim; ++i)
        os << ", " << rhs[i];
    return os << ")";
}

#endif
