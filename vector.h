#ifndef VECTOR_H
#define VECTOR_H

#include <array>
#include "Global.h"
#include "vector_like.h"

struct Vector : VectorLike<Vector, std::array<Long, dim>, dim> {
    using Container = std::array<Long, dim>;

    Container coords;

    Vector() { coords.fill(0); }

    Container &vector_container() { return coords; }
    const Container &vector_container() const { return coords; }
};

#endif
