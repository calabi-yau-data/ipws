#ifndef POINT_H
#define POINT_H

#include <array>
#include <gsl/gsl>
#include <utility>
#include <vector>
#include "global.h"
#include "vector_mixin.h"

struct Point : VectorMixin<Point, std::array<Long, dim>, dim> {
    using Container = std::array<Long, dim>;

    Container coords{};

    Container &vector_container() { return coords; }
    const Container &vector_container() const { return coords; }
};

// Returns a generating set of coordinate permutations that are symmetries of
// the given list of points. We could probably also consider symmetries that
// rely on changing the order of points, but this is not done here. Ihe indices
// in the pairs returned are in ascending orders.
std::vector<std::pair<unsigned, unsigned>> symmetries(
    gsl::span<const Point> points);

bool is_sorted(const Point &x,
               const std::vector<std::pair<unsigned, unsigned>> &checks);

#endif
