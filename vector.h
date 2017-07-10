#ifndef VECTOR_H
#define VECTOR_H

#include <array>
#include <cstddef>
#include <iostream>

#include "Global.h"
#include "vector_like.h"

struct VectorBase {
    std::array<Long, dim> coords;

    Long &vector_like_data(unsigned i) { return coords[i]; }
    const Long &vector_like_data(unsigned i) const { return coords[i]; }
};

struct Vector : VectorLike<VectorBase, Long, dim> {
    Vector() {}
    Vector(const Vector &v) { coords = v.coords; }
    Vector(const VectorBase &v) { coords = v.coords; }
};

#endif
