#ifndef WEIGHT_SYSTEM_BUILDER
#define WEIGHT_SYSTEM_BUILDER

#include <bitset>
#include <utility>
#include <vector>
#include "point.h"
#include "weight_system.h"

class WeightSystemBuilder {
public:
    WeightSystemBuilder();

    bool operator<(const WeightSystemBuilder &rhs) const;

    unsigned generator_count() const;

    WeightSystemBuilder restrict(const Point &x) const;
    bool average_if_nonzero(WeightSystem &q) const;
    bool allows(const Point &x) const;

    friend std::ostream &operator<<(std::ostream &os,
                                    const WeightSystemBuilder &rhs);

private:
    struct Generator {
        WeightSystem eq;
        std::bitset<dim> incidences;
    };
    std::vector<Generator>
        generators; // TODO: can there ever be more generators than dimensions?
    unsigned iteration;

    struct NoInit {
    };

    WeightSystemBuilder(NoInit) {}
};

bool leads_to_allowed_weightsystem(const Point &x);

#endif
