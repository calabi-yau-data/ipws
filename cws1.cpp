#include <iostream>
#include <set>
#include "config.h"
#include "stl_utils.h"
#include "stopwatch.h"
#include "vector.h"
#include "weight_system.h"
#include "weight_system_builder.h"

using std::array;
using std::cout;
using std::endl;
using std::set;
using std::setprecision;
using std::vector;

struct History {
    array<Vector, dim - 1> points;
    array<WeightSystem, dim> weight_systems;
    array<array<Long, dim - 1>, dim - 1> point_weight_system_distances;
};

struct Statistics {
    unsigned weight_systems_found;
};

__attribute__((noinline)) bool is_sorted(
    const Vector &x, const std::vector<std::pair<unsigned, unsigned>> &checks)
{
    for (const auto &check : checks)
        if (x.coords[check.first] < x.coords[check.second])
            return false;
    return true;
}

using WeightSystemCollection = set<WeightSystem>; // TODO: unordered_set?

bool add_maybe(WeightSystemCollection &weight_systems, WeightSystem ws)
{
    Long n = norm(ws);

    for (unsigned i = 0; i < dim; ++i) {
        if (!allow_weight_one_half &&
            2 * ws.weights[i] * r_numerator == n * r_denominator)
            return false;

        if (!allow_weight_one &&
            ws.weights[i] * r_numerator == n * r_denominator)
            return false;
    }

    if (!allow_weights_sum_one)
        for (unsigned i = 0; i < dim - 1; ++i)
            for (unsigned j = i + 1; j < dim; ++j)
                if ((ws.weights[i] + ws.weights[j]) * r_numerator ==
                    n * r_denominator)
                    return false;

    sort(ws);
    weight_systems.insert(ws);

    return true;
}

// TODO: verify this function
bool last_point_redundant2(const WeightSystemBuilder &builder, int n,
                           const History &history)
{
    auto &ws = history.weight_systems[n];

    auto points = WeightSystemPointsOn(ws);
    while (points.find_next()) {
        auto &x = points.get();

        if (!builder.allows(x))
            continue;

        // TODO: Do we have to check if x is linearly independent of the
        // other points?

        for (int i = 0; i < n; ++i) {
            Long diff = history.point_weight_system_distances[i][i] -
                        distance(history.weight_systems[i], x);
            if (diff > 0 ||
                (!disable_lex_compare && diff == 0 && x > history.points[i]))
                return true;
        }
    }

    return false;
}

// TODO: verify this function
bool last_point_redundant(int n, const History &history)
{
    Vector x = history.points[n];

    for (int i = 0; i < n; ++i) {
        Long rel = history.point_weight_system_distances[i][i] -
                   history.point_weight_system_distances[n][i];
        if (rel > 0)
            return true;
        if (!disable_lex_compare && rel == 0 && history.points[i] < x)
            return true;
    }

    for (int i = 0; i < n; ++i) {
        Vector x_other = history.points[i];
        Vector x_diff = x_other - x;

        Long v = gcd(x_diff);
        if (v != 1)
            x_diff /= v;

        bool all_positive = true;
        for (int j = 0; j < dim; ++j) {
            if (x_other.coords[j] + x_diff.coords[j] < 0) {
                all_positive = false;
                break;
            }
        }
        if (all_positive)
            return true;

        all_positive = true;
        for (int j = 0; j < dim; ++j) {
            if (x.coords[j] - x_diff.coords[j] < 0) {
                all_positive = false;
                break;
            }
        }
        if (all_positive)
            return true;
    }

    return false;
}

void rec(WeightSystemCollection &weight_systems,
         const WeightSystemBuilder &builder, int n, History &history,
         Statistics &statistics)
{
    WeightSystem ws{};
    if (!builder.average_if_nonzero(ws))
        return;

    history.weight_systems[n] = ws;

    if (n < dim - redundancy_check_skip_recursions &&
        last_point_redundant2(builder, n, history))
        return;

    if (add_maybe(weight_systems, ws))
        ++statistics.weight_systems_found;

    switch (n) {
    case dim - 2:
        // The following happens when redundancies are not checked
        assert(builder.generator_count() == 2);
        // TODO: special case for last iteration brings a big performance boost
        if (defer_last_recursion) {
            // builders.insert(builder.canonicalized());
            // ++buildersAddedCount;

            // if (buildersAddedCount % 10000 == 0)
            //     System.out.printf("%7.2f: builders: %d, unique: %d\n",
            //                       stopwatch.count(), buildersAddedCount,
            //                       builders.size());
            return;
        }
        break;
    case dim - 1:
        return;
    }

    auto symmetries = builder.symmetries();

    auto points = WeightSystemPointsBelow(ws);
    while (points.find_next()) {
        const Vector &x = points.get();

        if (!leads_to_allowed_weightsystem(x) ||
            (!debug_ignore_symmetries && !is_sorted(x, symmetries)))
            continue;

        history.points[n] = x;
        for (int i = 0; i < n + 1; ++i)
            history.point_weight_system_distances[n][i] =
                distance(history.weight_systems[i], x);

        if (last_point_redundant(n, history))
            continue;

        rec(weight_systems, builder.restrict(x), n + 1, history, statistics);
    }
}

int main()
{
    Stopwatch stopwatch{};

    History history{};
    WeightSystemCollection weight_systems{};
    Statistics statistics{};

    rec(weight_systems, WeightSystemBuilder{}, 0, history, statistics);

    if (defer_last_recursion) {
        // System.out.printf("%7.2f: builders: %d, unique: %d\n",
        //                   stopwatch.count(), buildersAddedCount,
        //                   builders.size());
        // System.out.printf("%7.2f: candidates: %d, unique: %d\n",
        //                   stopwatch.count(), weightSystemsAddedCount,
        //                   weightSystems.size());

        // for (WeightSystemBuilder builder : builders) {
        //     WeightSystem ws = builder.averageIfNonzero();

        //     WeightSystemBuilder.Symmetries symmetries = builder.symmetries();

        //     WeightSystem.PointsBelow points = ws.pointsBelow();
        //     while (points.findNext()) {
        //         Vector x = points.get();

        //         if (!WeightSystemBuilder.leadsToAllowedWeightsystem(x) ||
        //             (!debugIgnoreSymmetries && !symmetries.isSorted(x)))
        //             continue;

        //         WeightSystem ws2 = builder.restrict(x).averageIfNonzero();
        //         if (ws2 != null)
        //             addMaybe(ws2);
        //     }
        // }
    }

    cout << stopwatch
         << " - weight systems: " << statistics.weight_systems_found
         << ", unique: " << weight_systems.size() << endl;

    return 0;
}
