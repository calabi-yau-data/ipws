#ifndef WEIGHT_SYSTEM_H
#define WEIGHT_SYSTEM_H

#include <array>
#include "Global.h"
#include "vector.h"
#include "vector_like.h"

struct WeightSystem : VectorLike<WeightSystem, std::array<Long, dim>, dim> {
    using Container = std::array<Long, dim>;

    Container weights;

    WeightSystem() { weights.fill(0); }

    Container &vector_container() { return weights; }
    const Container &vector_container() const { return weights; }
};

class WeightSystemPointsBelow {
    WeightSystem q;
    Vector x;
    std::array<Long, dim> ax;

public:
    WeightSystemPointsBelow(const WeightSystem &q);
    const Vector &get() { return x; }
    bool find_next();
};

class WeightSystemPointsOn {
    WeightSystem q;
    Vector x;
    std::array<Long, dim> ax;

public:
    WeightSystemPointsOn(const WeightSystem &q);
    const Vector &get() { return x; }
    bool find_next();
};

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

#endif
