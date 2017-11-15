#ifndef WEIGHT_SYSTEM_H
#define WEIGHT_SYSTEM_H

#include <array>
#include "buffered_reader.h"
#include "buffered_writer.h"
#include "config.h"
#include "point.h"
#include "vector_mixin.h"

template <unsigned D>
struct WeightSystem : VectorMixin<WeightSystem<D>, std::array<Ring, D>, D> {
    using Container = std::array<Ring, D>;

    Container weights{};

    Container &vector_container() { return weights; }
    const Container &vector_container() const { return weights; }

    static constexpr size_t storage_size = D * sizeof(uint32_t);
};

namespace std {
template <unsigned D>
struct hash<WeightSystem<D>> {
    size_t operator()(WeightSystem<D> const &ws) const
    {
        size_t ret = 17280498655953172247u;
        for (const auto &w : ws.weights)
            ret = ret * 15004982843340118561u + std::hash<Ring>{}(w);
        return ret;
    }
};
}

class WeightSystemPointsBelow {
    WeightSystem<dim> q;
    Point x;
    std::array<Ring, dim> ax;

public:
    WeightSystemPointsBelow(const WeightSystem<dim> &q);
    const Point &get() { return x; }
    bool find_next();
};

class WeightSystemPointsOn {
    WeightSystem<dim> q;
    Point x;
    std::array<Ring, dim> ax;

public:
    WeightSystemPointsOn(const WeightSystem<dim> &q);
    const Point &get() { return x; }
    bool find_next();
};

Ring distance(const WeightSystem<dim> &ws, const Point &x);
bool good_weight_system(const WeightSystem<dim> &ws);

template <unsigned D>
void cancel(WeightSystem<D> &ws)
{
    if (D == 0)
        return;

    Ring gcd = std::abs(ws.weights[0]);

    for (unsigned i = 1; i < D; ++i)
        gcd = boost::math::gcd(gcd, ws.weights[i]);

    if (gcd != 1 && gcd != 0)
        ws /= gcd;
}

template <unsigned D>
void sort(WeightSystem<D> &ws)
{
    std::sort(ws.weights.begin(), ws.weights.end());
}

template <unsigned D>
Ring norm(const WeightSystem<D> &ws)
{
    return std::accumulate(ws.weights.begin(), ws.weights.end(), 0);
}

template <unsigned D>
void read(BufferedReader &f, WeightSystem<D> &ws)
{
    static_assert(WeightSystem<D>::storage_size == D * sizeof(int32_t),
                  "The constant WeightSystem::storage_size does not have the "
                  "right value");
    for (unsigned i = 0; i < D; ++i) {
        int32_t v;
        read(f, v);
        ws.weights[i] = v;
    }
}

template <unsigned D>
void write(BufferedWriter &f, const WeightSystem<D> &ws)
{
    for (unsigned i = 0; i < D; ++i) {
        auto v = ws.weights[i];
        assert(v >= 0 && v <= INT32_MAX);
        write32i(f, v);
    }
}

template <unsigned D>
void read_varint(BufferedReader &f, WeightSystem<D> &ws)
{
    for (unsigned i = 0; i < D; ++i) {
        auto v = read_varint(f);
        assert(v <= std::numeric_limits<Ring>::max());
        ws.weights[i] = static_cast<Ring>(v);
    }
}

template <unsigned D>
void write_varint(BufferedWriter &f, const WeightSystem<D> &ws)
{
    for (unsigned i = 0; i < D; ++i)
        write_varint(f, ws.weights[i]);
}

// Returns the weight system q that is a linear combination of q1 and q2
// such that its distance to x is zero. Equivalently:
// 1) distance(q, x) == 0
// 2) for any y, distance(q1, y) == 0 && distance(q2, y) == 0 implies
//    distance(q, y) == 0
const WeightSystem<dim> intersect(const WeightSystem<dim> &q1,
                                  const WeightSystem<dim> &q2, const Point &x);

bool has_ip(const WeightSystem<dim> &ws);

#endif
