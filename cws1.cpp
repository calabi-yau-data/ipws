#include <iostream>
#include <set>

#include "config.h"
#include "stl_utils.h"
#include "stopwatch.h"
#include "weight_system.h"
#include "weight_system_builder.h"
#include "vector.h"

FILE *inFILE;
FILE *outFILE;

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
    Long norm = std::accumulate(ws.a.begin(), ws.a.end(), 0);

    for (unsigned i = 0; i < dim; ++i) {
        if (!allow_weight_one_half && 2 * ws.a[i] * r_numerator == norm * r_denominator)
            return;

        if (!allow_weight_one && ws.a[i] * r_numerator == norm * r_denominator)
            return;
    }

    if (!allow_weights_sum_one)
        for (unsigned i = 0; i < dim - 1; ++i)
            for (unsigned j = i + 1; j < dim; ++j)
                if ((ws.a[i] + ws.a[j]) * r_numerator == norm * r_denominator)
                    return;

    ++count;

    ws.sort();
    weight_systems.insert(ws);
}

void rec(const WeightSystemBuilder &builder) {
    WeightSystem ws{};
    if (!builder.average_if_nonzero(ws))
        return;

    // // it is not economical to do this on the last two recursion levels
    // if (builder.iteration() < dim - 2) {
    //     bool skip = false;

    //     auto points = points_on(ws);
    //     while (points.find_next()) {
    //         const Vector &x = points.get();

    //         if (!builder.allows(x))
    //             continue;

    //         // TODO: Do we have to check if x is linearly independent of the
    //         // other points?
    //         for (unsigned i = 0; i < n; ++i) {
    //             Long diff = X.x_inner_q[i + 1][i] - x * X.qs[i];
    //             if (diff > 0 || (diff == 0 && lex_cmp(x, X.x[i + 1]) > 0)) {
    //                 skip = true;
    //                 return;
    //             }
    //         }
    //     });

    //     if (skip)
    //         return;
    // }

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

        if (!WeightSystemBuilder::leads_to_allowed_weightsystem(x, r_numerator,
                                                                r_denominator) ||
            (!debug_ignore_symmetries && !is_sorted(x, symmetries)))
            continue;

        rec(builder.restrict(x));
    }
}

int main() {
    Stopwatch stopwatch{};

    rec(WeightSystemBuilder{});

    cout << stopwatch.count() << ": "
         << weight_systems.size() << "/" << count << endl;

    return 0;
}

// TODO: Acknowledge Roman
