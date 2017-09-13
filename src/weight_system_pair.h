#ifndef WEIGHT_SYSTEM_PAIR_H
#define WEIGHT_SYSTEM_PAIR_H

#include <algorithm>
#include <gsl/gsl>
#include <utility>
#include <vector>
#include "weight_system.h"

using WeightSystemPair = std::pair<WeightSystem<dim>, WeightSystem<dim>>;

namespace std {
template <>
struct hash<WeightSystemPair> {
    size_t operator()(WeightSystemPair const &pair) const
    {
        size_t ret = 13968802923818332081u;
        ret = ret * 14053336887773799961u +
              std::hash<WeightSystem<dim>>{}(pair.first);
        ret = ret * 14053336887773799961u +
              std::hash<WeightSystem<dim>>{}(pair.second);
        return ret;
    }
};
}

const WeightSystemPair canonicalize(const WeightSystemPair &pair);
const WeightSystem<dim> average(const WeightSystemPair &pair);
bool restrict(const WeightSystemPair &pair, const Point &x,
              WeightSystem<dim> &ws);

// Returns a generating set of coordinate permutations that are symmetries of
// the pair of weight systems, while allowing exchange of weight systems. The
// indices in the pairs returned are in ascending orders.
// std::vector<std::pair<unsigned, unsigned>> symmetries(
//     const WeightSystemPair &pair);

#endif
