#ifndef WEIGHT_SYSTEM_BUILDER
#define WEIGHT_SYSTEM_BUILDER

#include "point.h"
#include "weight_system.h"
#include "weight_system_pair.h"
#include <bitset>
#include <utility>
#include <vector>

class WeightSystemBuilder {
public:
    WeightSystemBuilder();

    bool operator<(const WeightSystemBuilder &rhs) const;

    unsigned generator_count() const;

    const WeightSystemBuilder restrict(const Point &x) const;
    bool average_if_nonzero(WeightSystem<dim> &q) const;
    bool allows(const Point &x) const;

    friend std::ostream &operator<<(std::ostream &os,
                                    const WeightSystemBuilder &rhs);

    const WeightSystemPair to_pair() const;

private:
    struct Generator {
        WeightSystem<dim> eq;
        std::bitset<dim> incidences;
    };
    std::vector<Generator> generators;
    unsigned iteration;

    struct NoInit {};

    WeightSystemBuilder(NoInit) {}
};

bool leads_to_allowed_weightsystem(const Point &x);

#endif
