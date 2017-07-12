#ifndef VECTOR_H
#define VECTOR_H

#include <array>
#include "global.h"
#include "vector_mixin.h"

struct Point : VectorMixin<Point, std::array<Long, dim>, dim> {
    using Container = std::array<Long, dim>;

    Container coords{};

    Container &vector_container() { return coords; }
    const Container &vector_container() const { return coords; }
};

#endif
