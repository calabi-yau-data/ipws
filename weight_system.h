#ifndef WEIGHT_SYSTEM_H
#define WEIGHT_SYSTEM_H

#include <array>
#include "global2.h"
#include "point.h"
#include "vector_mixin.h"

struct WeightSystem : VectorMixin<WeightSystem, std::array<Long, dim>, dim> {
    using Container = std::array<Long, dim>;

    Container weights{};

    Container &vector_container() { return weights; }
    const Container &vector_container() const { return weights; }
};

class WeightSystemPointsBelow {
    WeightSystem q;
    Point x;
    std::array<Long, dim> ax;

public:
    WeightSystemPointsBelow(const WeightSystem &q);
    const Point &get() { return x; }
    bool find_next();
};

class WeightSystemPointsOn {
    WeightSystem q;
    Point x;
    std::array<Long, dim> ax;

public:
    WeightSystemPointsOn(const WeightSystem &q);
    const Point &get() { return x; }
    bool find_next();
};

Long distance(const WeightSystem &ws, const Point &x);
void cancel(WeightSystem &ws);
void sort(WeightSystem &ws);
Long norm(const WeightSystem &ws);

// Returns the weight system q that is a linear combination of q1 and q2
// such that its distance to x is zero. Equivalently:
// 1) q.distance_from(x) == 0
// 2) q1.distance_from(y) == 0 && q2.distance_from(y) == 0 implies
//    q.distance_from(y) == 0 for all y
WeightSystem intersect(const WeightSystem &q1, const WeightSystem &q2,
                       const Point &x);

#endif
