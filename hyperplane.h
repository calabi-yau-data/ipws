#ifndef HYPERPLANE_H
#define HYPERPLANE_H

#include "Global.h"
#include "vector_old.h"

#include <experimental/numeric>

// Represents the hyperplane given by all x such that x*a + c = 0.
struct Hyperplane {
    Vector a;
    Long c;

    Hyperplane() {}
    Hyperplane(const Vector &a, Long c) : a(a), c(c) {}

    bool operator!=(const Hyperplane &rhs) const
    {
        if (c != rhs.c)
            return true;

        for (size_t i = 0; i < dim; ++i)
            if (a[i] != rhs.a[i])
                return true;

        return false;
    }

    bool operator<(const Hyperplane &rhs) const
    {
        if (c != rhs.c)
            return c < rhs.c;

        for (size_t i = 0; i < dim; ++i)
            if (a[i] != rhs.a[i])
                return a[i] < rhs.a[i];

        return false;
    }

    Hyperplane operator-() const { return Hyperplane(-a, -c); }
    Long operator*(const Vector &rhs) const
    {
        Long ret{c};
        for (size_t i = 0; i < dim; ++i)
            ret += rhs[i] * a[i];
        return ret;
    }

    friend Long operator*(const Vector &lhs, const Hyperplane &rhs)
    {
        Long ret{rhs.c};
        for (size_t i = 0; i < dim; ++i)
            ret += lhs[i] * rhs.a[i];
        return ret;
    }

    friend std::ostream &operator<<(std::ostream &os, const Hyperplane &rhs)
    {
        os << "(" << rhs.c << ";";
        if (dim != 0)
            os << " " << rhs.a[0];
        for (size_t i = 1; i < dim; ++i)
            os << ", " << rhs.a[i];
        return os << ")";
    }

    void cancel()
    {
        Long gcd = std::abs(c);

        for (size_t i = 0; i < dim; ++i)
            gcd = std::experimental::gcd(gcd, a[i]);

        if (gcd != 1 && gcd != 0) {
            a /= gcd;
            c /= gcd;
        }
    }

    // Returns the hyperplane through v and the intersection of q1 and q2.
    friend Hyperplane intersect(const Hyperplane &q1, const Hyperplane &q2,
                                const Vector &v)
    {
        Long e1 = v * q1;
        Long e2 = v * q2;

        Long gcd = std::experimental::gcd(e1, e2);
        e1 /= gcd;
        e2 /= gcd;

        Hyperplane ret{};
        ret.a = e2 * q1.a - e1 * q2.a;
        ret.c = e2 * q1.c - e1 * q2.c;

        ret.cancel();
        return ret;
    }
};

#endif
