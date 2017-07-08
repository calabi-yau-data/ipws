#ifndef WEIGHT_SYSTEM_H
#define WEIGHT_SYSTEM_H

#include "Global.h"
#include "vector.h"

#include <algorithm>
#include <experimental/numeric>
#include <numeric>

#include "config.h"

class WeightSystemPointsBelow;

struct WeightSystem {
    Vector weights;

    int compare(const WeightSystem &rhs) const {
        for (unsigned i = 0; i < dim; ++i) {
            if (weights[i] < rhs.weights[i])
                return -1;
            if (weights[i] > rhs.weights[i])
                return 1;
        }

        return 0;
    }

    bool operator==(const WeightSystem &rhs) const { return compare(rhs) == 0; }
    bool operator!=(const WeightSystem &rhs) const { return compare(rhs) != 0; }
    bool operator<(const WeightSystem &rhs) const { return compare(rhs) < 0; }
    bool operator>(const WeightSystem &rhs) const { return compare(rhs) > 0; }

    WeightSystem operator-() const { return WeightSystem{-weights}; }

    Long distance_to(const Vector &v) const {
        Long ret = 0;
        for (unsigned i = 0; i < dim; ++i)
            ret += (v[i] * r_numerator - r_denominator) * weights[i];
        return ret;
    }

    friend std::ostream &operator<<(std::ostream &os, const WeightSystem &rhs) {
        os << "(";
        if (dim != 0)
            os << rhs.weights[0];
        for (unsigned i = 1; i < dim; ++i)
            os << ", " << rhs.weights[i];
        return os << ")";
    }

    void cancel() {
        if (dim == 0)
            return;

        Long gcd = std::abs(weights[0]);

        for (unsigned i = 1; i < dim; ++i)
            gcd = std::experimental::gcd(gcd, weights[i]);

        if (gcd != 1 && gcd != 0)
            weights /= gcd;
    }

    // Returns the hyperplane through v and the intersection of q1 and q2.
    friend WeightSystem intersect(const WeightSystem &q1,
                                  const WeightSystem &q2, const Vector &v) {
        Long e1 = q1.distance_to(v);
        Long e2 = q2.distance_to(v);

        Long gcd = std::experimental::gcd(e1, e2);
        e1 /= gcd;
        e2 /= gcd;

        WeightSystem ret{};
        if (e1 > 0)
            ret.weights = e1 * q2.weights - e2 * q1.weights;
        else
            ret.weights = e2 * q1.weights - e1 * q2.weights;

        ret.cancel();
        return ret;
    }

    void sort() {
        std::sort(weights.begin(), weights.end());
    }

    Long norm() const {
        return std::accumulate(weights.begin(), weights.end(), 0);
    }
};

class WeightSystemPointsBelow {
    WeightSystem q;
    Vector x;
    std::array<Long, dim> ax;
public:
    WeightSystemPointsBelow(const WeightSystem &q) : q{q} {
        x.fill(0);
        ax.fill(-q.norm() * r_denominator);

        x[dim - 1] -= 1;
        ax[dim - 1] -= q.weights[dim - 1] * r_numerator;
    }

    const Vector &get() {
        return x;
    }

    __attribute__ ((noinline))
    bool find_next() {
        int k = dim - 1;
        while (ax[k] + q.weights[k] * r_numerator >= 0) {
            if (k == 0)
                return false;
            x[k] = 0;
            --k;
        }

        x[k]++;
        ax[k] += q.weights[k] * r_numerator;
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
        assert(false); // TODO: not yet shifted

        std::fill_n(x.begin(), dim - 1, 0);
        std::fill_n(ax.begin(), dim - 1, 0);

        x[dim - 1] = -1;
        ax[dim - 1] = -q.weights[dim - 1];
    }

    const Vector &get() {
        return x;
    }

    // this is just a quick adaptation and can probably be optimized
    __attribute__ ((noinline))
    bool find_next() {
        while (true) {
            int k = dim - 1;
            while (ax[k] + q.weights[k] > 0) {
                if (k == 0)
                    return false;
                x[k] = 0;
                --k;
            }

            x[k]++;
            ax[k] += q.weights[k];
            for (int i = k + 1; i < dim; ++i)
                ax[i] = ax[k];

            if (ax[k] == 0)
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
