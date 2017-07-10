#ifndef WEIGHT_SYSTEM_BUILDER
#define WEIGHT_SYSTEM_BUILDER

#include <algorithm>
#include <bitset>
#include <experimental/numeric>
#include <vector>
#include "stl_utils.h"
#include "vector.h"
#include "weight_system.h"

class WeightSystemBuilder {
public:
    WeightSystemBuilder();

    bool operator<(const WeightSystemBuilder &rhs) const;

    unsigned generator_count() const;

    // the indices in the pairs are in ascending orders
    std::vector<std::pair<unsigned, unsigned>> symmetries() const;
    bool has_symmetry(unsigned idx1, unsigned idx2) const;

    WeightSystemBuilder restrict(const Vector &x) const;
    bool average_if_nonzero(WeightSystem &q) const;
    bool allows(const Vector &x) const;

    static bool leads_to_allowed_weightsystem(const Vector &x, Long r_numerator,
                                              Long r_denominator);

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

#endif
