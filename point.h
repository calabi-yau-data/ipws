#ifndef POINT_H
#define POINT_H

#include <array>
#include <vector>
#include "global.h"
#include "vector_mixin.h"

struct Point : VectorMixin<Point, std::array<Long, dim>, dim> {
    using Container = std::array<Long, dim>;

    Container coords{};

    Container &vector_container() { return coords; }
    const Container &vector_container() const { return coords; }
};

bool is_sorted(const Point &x,
               const std::vector<std::pair<unsigned, unsigned>> &checks);

#endif
