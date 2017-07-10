#ifndef VECTOR_LIKE_H
#define VECTOR_LIKE_H

#include <iostream>

template <class T, class Container, unsigned D>
class VectorLike {
    using R = typename Container::value_type;

public:
    friend int compare(const T &a, const T &b)
    {
        for (unsigned i = 0; i < D; ++i) {
            if (a.vector_container()[i] < b.vector_container()[i])
                return -1;
            if (a.vector_container()[i] > b.vector_container()[i])
                return 1;
        }
        return 0;
    }

    friend bool operator==(const T &lhs, const T &rhs)
    {
        return compare(lhs, rhs) == 0;
    }
    friend bool operator!=(const T &lhs, const T &rhs)
    {
        return compare(lhs, rhs) != 0;
    }
    friend bool operator<(const T &lhs, const T &rhs)
    {
        return compare(lhs, rhs) < 0;
    }
    friend bool operator>(const T &lhs, const T &rhs)
    {
        return compare(lhs, rhs) > 0;
    }

    friend T &operator/=(T &lhs, const R &rhs)
    {
        for (unsigned i = 0; i < D; ++i)
            lhs.vector_container()[i] *= rhs;
        return lhs;
    }

    friend T operator-(const T &lhs, const T &rhs)
    {
        T ret{};
        for (unsigned i = 0; i < D; ++i)
            ret.vector_container()[i] =
                lhs.vector_container()[i] - rhs.vector_container()[i];
        return ret;
    }

    friend T operator*(const R &lhs, const T &rhs)
    {
        T ret{};
        for (unsigned i = 0; i < D; ++i)
            ret.vector_container()[i] = rhs.vector_container()[i] * lhs;
        return ret;
    }

    friend std::ostream &operator<<(std::ostream &os, const T &rhs)
    {
        os << "(";
        if (D != 0)
            os << rhs.vector_container()[0];
        for (unsigned i = 1; i < D; ++i)
            os << ", " << rhs.vector_container()[i];
        return os << ")";
    }
};

#endif
