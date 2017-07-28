#ifndef WEIGHT_SYSTEM_H
#define WEIGHT_SYSTEM_H

#include <array>
#include "buffered_reader.h"
#include "buffered_writer.h"
#include "config.h"
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
        size_t ret = 17280498655953172247u;
        for (const auto &w : ws.weights)
            ret = ret * 15004982843340118561u + std::hash<Ring>{}(w);
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
void read(BufferedReader &f, WeightSystem &ws);
void write(BufferedWriter &f, const WeightSystem &ws);
void read_varint(BufferedReader &f, WeightSystem &ws);
void write_varint(BufferedWriter &f, const WeightSystem &ws);

// Returns the weight system q that is a linear combination of q1 and q2
// such that its distance to x is zero. Equivalently:
// 1) distance(q, x) == 0
// 2) for any y, distance(q1, y) == 0 && distance(q2, y) == 0 implies
//    distance(q, y) == 0
const WeightSystem intersect(const WeightSystem &q1, const WeightSystem &q2,
                             const Point &x);

bool has_ip(const WeightSystem &ws);

constexpr unsigned weight_system_storage_size = dim * sizeof(uint32_t);

#endif
