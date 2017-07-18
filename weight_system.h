#ifndef WEIGHT_SYSTEM_H
#define WEIGHT_SYSTEM_H

#include <array>
#include "config.h"
#include "file.h"
#include "point.h"
#include "vector_mixin.h"

struct WeightSystem : VectorMixin<WeightSystem, std::array<Ring, dim>, dim> {
    using Container = std::array<Ring, dim>;

    Container weights{};

    Container &vector_container() { return weights; }
    const Container &vector_container() const { return weights; }
};

namespace std {
template <>
struct hash<WeightSystem> {
    size_t operator()(WeightSystem const &ws) const
    {
        size_t ret = 237037;
        for (const auto &w : ws.weights)
            ret = ret * 92647 + std::hash<Ring>{}(w);
        return ret;
    }
};
}

class WeightSystemPointsBelow {
    WeightSystem q;
    Point x;
    std::array<Ring, dim> ax;

public:
    WeightSystemPointsBelow(const WeightSystem &q);
    const Point &get() { return x; }
    bool find_next();
};

class WeightSystemPointsOn {
    WeightSystem q;
    Point x;
    std::array<Ring, dim> ax;

public:
    WeightSystemPointsOn(const WeightSystem &q);
    const Point &get() { return x; }
    bool find_next();
};

Ring distance(const WeightSystem &ws, const Point &x);
void cancel(WeightSystem &ws);
void sort(WeightSystem &ws);
Ring norm(const WeightSystem &ws);
bool good_weight_system(const WeightSystem &ws);
void read(File &f, WeightSystem &ws);
void write(File &f, const WeightSystem &ws);

// Returns the weight system q that is a linear combination of q1 and q2
// such that its distance to x is zero. Equivalently:
// 1) q.distance_from(x) == 0
// 2) q1.distance_from(y) == 0 && q2.distance_from(y) == 0 implies
//    q.distance_from(y) == 0 for all y
const WeightSystem intersect(const WeightSystem &q1, const WeightSystem &q2,
                             const Point &x);

bool has_ip(const WeightSystem &ws);

#endif
