#ifndef WEIGHT_SYSTEM_PAIR_H
#define WEIGHT_SYSTEM_PAIR_H

#include <algorithm>
#include <utility>
#include "weight_system.h"

using WeightSystemPair = std::pair<WeightSystem, WeightSystem>;

WeightSystemPair canonicalize(const WeightSystemPair &pair);

#endif
