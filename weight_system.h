#ifndef HYPERPLANE_H
#define HYPERPLANE_H

#include "Global.h"
#include "vector.h"

#include <algorithm>
#include <experimental/numeric>

#include "config.h"

class WeightSystemPointsBelow;

// Represents the weight system q = a / c
struct WeightSystem {
    Vector a;
    Long c;

    int compare(const WeightSystem &rhs) const {
        if (c < rhs.c)
            return -1;
        if (c > rhs.c)
            return 1;

        for (size_t i = 0; i < dim; ++i) {
            if (a[i] < rhs.a[i])
                return -1;
            if (a[i] > rhs.a[i])
                return 1;
        }

        return 0;
    }

    bool operator==(const WeightSystem &rhs) const { return compare(rhs) == 0; }
    bool operator!=(const WeightSystem &rhs) const { return compare(rhs) != 0; }
    bool operator<(const WeightSystem &rhs) const { return compare(rhs) < 0; }
    bool operator>(const WeightSystem &rhs) const { return compare(rhs) > 0; }

    WeightSystem operator-() const { return WeightSystem{-a, -c}; }

    Long apply_to(const Vector &v) const {
        Long ret{-c};
        for (size_t i = 0; i < dim; ++i)
            ret += v[i] * a[i];
        return ret;
    }

    friend std::ostream &operator<<(std::ostream &os, const WeightSystem &rhs) {
        os << "(" << rhs.c << ";";
        if (dim != 0)
            os << " " << rhs.a[0];
        for (size_t i = 1; i < dim; ++i)
            os << ", " << rhs.a[i];
        return os << ")";
    }

    void cancel() {
        Long gcd = 1;
        if (c > 0)
            gcd = c;
        else if (c < 0)
            gcd = -c;

        for (size_t i = 0; i < dim; ++i)
            gcd = std::experimental::gcd(gcd, a[i]);

        if (gcd != 1) {
            a /= gcd;
            c /= gcd;
        }
    }

    // Returns the hyperplane through v and the intersection of q1 and q2.
    friend WeightSystem intersect(const WeightSystem &q1,
                                  const WeightSystem &q2, const Vector &v) {
        Long e1 = q1.apply_to(v);
        Long e2 = q2.apply_to(v);

        Long gcd = std::experimental::gcd(e1, e2);
        e1 /= gcd;
        e2 /= gcd;

        WeightSystem ret{};
        ret.a = e2 * q1.a - e1 * q2.a;
        ret.c = e2 * q1.c - e1 * q2.c;

        ret.cancel();
        return ret;
    }

    void sort() {
        std::sort(a.begin(), a.end());
    }
};

class WeightSystemPointsBelow {
    WeightSystem q;
    Vector x;
    std::array<Long, dim> ax;
public:
    WeightSystemPointsBelow(const WeightSystem &q) : q{q} {
        std::fill_n(x.begin(), dim - 1, 0);
        std::fill_n(ax.begin(), dim - 1, 0);

        x[dim - 1] = -1;
        ax[dim - 1] = -q.a[dim - 1];
    }

    const Vector &get() {
        return x;
    }

    __attribute__ ((noinline))
    bool find_next() {
        int k = dim - 1;
        while (ax[k] + q.a[k] >= q.c) {
            if (k == 0)
                return false;
            x[k] = 0;
            --k;
        }

        x[k]++;
        ax[k] += q.a[k];
        for (int i = k + 1; i < dim; ++i)
            ax[i] = ax[k];

        return true;
    }
};

class WeightSystemPointsOn {
    WeightSystem q;
    Vector x;
    std::array<Long, dim> ax;
public:
    WeightSystemPointsOn(const WeightSystem &q) : q{q} {
        std::fill_n(x.begin(), dim - 1, 0);
        std::fill_n(ax.begin(), dim - 1, 0);

        x[dim - 1] = -1;
        ax[dim - 1] = -q.a[dim - 1];
    }

    const Vector &get() {
        return x;
    }

    // this is just a quick adaptation and can probably be optimized
    __attribute__ ((noinline))
    bool find_next() {
        while (true) {
            int k = dim - 1;
            while (ax[k] + q.a[k] > q.c) {
                if (k == 0)
                    return false;
                x[k] = 0;
                --k;
            }

            x[k]++;
            ax[k] += q.a[k];
            for (int i = k + 1; i < dim; ++i)
                ax[i] = ax[k];

            if (ax[k] == q.c)
                return true;
        };
    }
};

__attribute__ ((noinline))
WeightSystemPointsBelow points_below(WeightSystem &q) {
    return WeightSystemPointsBelow{q};
}

__attribute__ ((noinline))
WeightSystemPointsOn points_on(WeightSystem &q) {
    return WeightSystemPointsOn{q};
}

#endif
