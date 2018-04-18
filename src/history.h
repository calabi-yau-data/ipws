#ifndef HISTORY_H
#define HISTORY_H

#include <array>
#include "point.h"
#include "weight_system.h"
#include "weight_system_builder.h"

struct History {
    template <class T, size_t D>
    using array = std::array<T, D>;

    array<Point, dim - 1> points{};
    array<WeightSystem<dim>, dim> weight_systems{};
    array<array<Ring, dim - 1>, dim - 1> point_weight_system_distances{};
};

bool last_point_redundant_quick(int n, const History &history);
bool last_point_redundant(const WeightSystemBuilder &builder, int n,
                          const History &history);

#endif
