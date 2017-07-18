#ifndef WEIGHT_SYSTEM_PAIR_H
#define WEIGHT_SYSTEM_PAIR_H

#include <algorithm>
#include <gsl/gsl>
#include <utility>
#include <vector>
#include "weight_system.h"

using WeightSystemPair = std::pair<WeightSystem, WeightSystem>;

const WeightSystemPair canonicalize(const WeightSystemPair &pair);
const WeightSystem average(const WeightSystemPair &pair);
bool restrict(const WeightSystemPair &pair, const Point &x, WeightSystem &ws);

// Returns a generating set of coordinate permutations that are symmetries of
// the pair of weight systems, while allowing exchange of weight systems. The
// indices in the pairs returned are in ascending orders.
// std::vector<std::pair<unsigned, unsigned>> symmetries(
//     const WeightSystemPair &pair);

#endif
