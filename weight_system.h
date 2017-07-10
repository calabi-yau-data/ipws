#ifndef WEIGHT_SYSTEM_H
#define WEIGHT_SYSTEM_H

#include "Global.h"
#include "vector.h"
#include "vector_like.h"

#include "config.h"

struct WeightSystem;

Long distance(const WeightSystem &ws, const Vector &x);
void cancel(WeightSystem &ws);
void sort(WeightSystem &ws);
Long norm(const WeightSystem &ws);

// Returns the weight system q that is a linear combination of q1 and q2
// such that its distance to x is zero. Equivalently:
// 1) q.distance_from(x) == 0
// 2) q1.distance_from(y) == 0 && q2.distance_from(y) == 0 implies
//    q.distance_from(y) == 0 for all y
WeightSystem intersect(const WeightSystem &q1, const WeightSystem &q2,
                       const Vector &x);

struct WeightSystem : VectorLike<WeightSystem, std::array<Long, dim>, dim> {
    using Container = std::array<Long, dim>;

    Container weights;

    Container &vector_container() { return weights; }
    const Container &vector_container() const { return weights; }
};

class WeightSystemPointsBelow {
    WeightSystem q;
    Vector x;
    std::array<Long, dim> ax;

public:
    WeightSystemPointsBelow(const WeightSystem &q) : q{q}
    {
        x.coords.fill(0);
        ax.fill(-norm(q) * r_denominator);

        x.coords[dim - 1] -= 1;
        ax[dim - 1] -= q.weights[dim - 1] * r_numerator;
    }

    const Vector &get() { return x; }
    bool find_next()
    {
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
    WeightSystemPointsOn(const WeightSystem &q) : q{q}
    {
        x.coords.fill(0);
        ax.fill(-norm(q) * r_denominator);

        x.coords[dim - 1] -= 1;
        ax[dim - 1] -= q.weights[dim - 1] * r_numerator;
    }

    const Vector &get() { return x; }
    bool find_next()
    {
        while (true) {
            int k = dim - 1;
            while (ax[k] + q.weights[k] * r_numerator > 0) {
                if (k == 0)
                    return false;
                x.coords[k] = 0;
                --k;
            }

            x.coords[k]++;
            ax[k] += q.weights[k] * r_numerator;
            for (int i = k + 1; i < dim; ++i)
                ax[i] = ax[k];

            if (ax[k] == 0)
                return true;
        }
    }
};

#endif
