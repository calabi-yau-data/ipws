#include <iostream>
#include <set>

#include "stl_utils.h"
#include "stopwatch.h"
#include "weight_system.h"
#include "weight_system_builder.h"
#include "vector.h"

FILE *inFILE;
FILE *outFILE;

const Long r_numerator = 3;
const Long r_denominator = 2;

using std::array;
using std::cout;
using std::endl;
using std::set;
using std::vector;

__attribute__ ((noinline))
bool is_sorted(const Vector &x, const std::vector<std::pair<unsigned, unsigned>> &checks) {
    for (const auto &check : checks)
        if (x[check.first] < x[check.second])
            return false;
    return true;
}

unsigned count = 0;
set<WeightSystem> weight_systems{};

void add_maybe(WeightSystem ws) {
    for (unsigned i = 0; i < dim; ++i) {
        // Skip weight systems containing a weight of 1/2
        if (2 * ws.a[i] == ws.c)
            return;

        // Skip weight systems containing a weight of 1
        if (ws.a[i] == ws.c)
            return;
    }

    // Skip weight systems containing two weights with a sum of 1
    for (unsigned i = 0; i < dim - 1; ++i)
        for (unsigned j = i + 1; j < dim; ++j)
            if (ws.a[i] + ws.a[j] == ws.c)
                return;

    ++count;

    ws.sort();
    weight_systems.insert(ws);
}

void rec(const WeightSystemBuilder &builder) {
    WeightSystem ws{};
    if (!builder.average_if_nonzero(ws))
        return;

    add_maybe(ws);

    switch (builder.iteration()) {
    case dim - 2:
        assert(builder.generator_count() == 2);
        // TODO: special case for last iteration brings a big performance boost
        break;
    case dim - 1:
        return; // TODO: is this not needed in the old code?
    }

    auto symmetries = builder.symmetries();

    auto points = points_below(ws);
    while (points.find_next()) {
        const Vector &x = points.get();

        if (!WeightSystemBuilder::gives_good_weightsystem(x, r_numerator,
                                                          r_denominator) ||
            !is_sorted(x, symmetries))
            continue;

        rec(builder.restrict(x));
    }
}

int main() {
    Stopwatch stopwatch{};

    rec(WeightSystemBuilder{r_numerator, r_denominator});

    cout << stopwatch.count() << ": "
         << weight_systems.size() << "/" << count << endl;

    return 0;
}

// TODO: Acknowledge Roman
