#ifndef WEIGHT_SYSTEM_H
#define WEIGHT_SYSTEM_H

#include "Global.h"
#include "vector.h"
#include "vector_like.h"

#include <algorithm>
#include <experimental/numeric>
#include <numeric>

#include "config.h"

class WeightSystemPointsBelow;

struct WeightSystemBase {
    std::array<Long, dim> weights;

    Long &vector_like_data(unsigned i) { return weights[i]; }
    const Long &vector_like_data(unsigned i) const { return weights[i]; }
};

struct WeightSystem : VectorLike<WeightSystemBase, Long, dim> {
    WeightSystem() {}
    WeightSystem(const WeightSystem &q) { weights = q.weights; }
    WeightSystem(const WeightSystemBase &q) { weights = q.weights; }

    Long distance_from(const Vector &v) const {
        Long ret = 0;
        for (unsigned i = 0; i < dim; ++i)
            ret += (v.coords[i] * r_numerator - r_denominator) * weights[i];
        return ret;
    }

    void cancel() {
        if (dim == 0)
            return;

        Long gcd = std::abs(weights[0]);

        for (unsigned i = 1; i < dim; ++i)
            gcd = std::experimental::gcd(gcd, weights[i]);

        if (gcd != 1 && gcd != 0)
            *this /= gcd;
    }

    // Returns the hyperplane through v and the intersection of q1 and q2.
    friend WeightSystem intersect(const WeightSystem &q1,
                                  const WeightSystem &q2, const Vector &v) {
        Long e1 = q1.distance_from(v);
        Long e2 = q2.distance_from(v);

        Long gcd = std::experimental::gcd(e1, e2);
        e1 /= gcd;
        e2 /= gcd;

        WeightSystem ret{};
        if (e1 < 0)
            ret = e1 * q2 - e2 * q1;
        else
            ret = e2 * q1 - e1 * q2;

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
        x.coords.fill(0);
        ax.fill(-q.norm() * r_denominator);

        x.coords[dim - 1] -= 1;
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
            x.coords[k] = 0;
            --k;
        }

        x.coords[k]++;
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

        std::fill_n(x.coords.begin(), dim - 1, 0);
        std::fill_n(ax.begin(), dim - 1, 0);

        x.coords[dim - 1] = -1;
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
                x.coords[k] = 0;
                --k;
            }

            x.coords[k]++;
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
