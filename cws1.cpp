#include <gsl/gsl>
#include <iostream>
#include <set>
#include "config.h"
#include "point.h"
#include "stl_utils.h"
#include "stopwatch.h"
#include "weight_system.h"
#include "weight_system_builder.h"
#include "weight_system_pair.h"

using gsl::span;
using std::array;
using std::cout;
using std::endl;
using std::pair;
using std::set;
using std::setprecision;
using std::vector;

set<WeightSystemPair> pairs; // TODO: unordered_set? TODO: not global

struct History {
    array<Point, dim - 1> points;
    array<WeightSystem, dim> weight_systems;
    array<array<Long, dim - 1>, dim - 1> point_weight_system_distances;
};

struct Statistics {
    unsigned weight_systems_found;
    unsigned final_pairs_found;
};

__attribute__((noinline)) bool is_sorted(
    const Point &x, const std::vector<std::pair<unsigned, unsigned>> &checks)
{
    for (const auto &check : checks)
        if (x.coords[check.first] < x.coords[check.second])
            return false;
    return true;
}

using WeightSystemCollection = set<WeightSystem>; // TODO: unordered_set?

bool points_have_symmetry(span<Point> points, unsigned idx1, unsigned idx2)
{
    // We could probably also consider symmetries that rely on changing the
    // order of points, but this is not done here.

    for (const auto &p : points)
        if (p.coords[idx1] != p.coords[idx2])
            return false;
    return true;
}

// the indices in the pairs are in ascending orders
vector<pair<unsigned, unsigned>> points_symmetries(span<Point> points)
{
    vector<pair<unsigned, unsigned>> ret{};

    array<bool, dim> done{};
    done.fill(false);

    for (unsigned i = 0; i < dim - 1; ++i) {
        if (done[i])
            continue;

        unsigned k = i;
        for (unsigned j = i + 1; j < dim; ++j) {
            if (done[j])
                continue;

            if (points_have_symmetry(points, i, j)) {
                ret.push_back(pair<unsigned, unsigned>(k, j));
                k = j;
                done[j] = true;
            }
        }
    }

    return ret;
}

bool good_weight_system(const WeightSystem &ws)
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

    return true;
}

void add_maybe(WeightSystemCollection &weight_systems, WeightSystem ws,
               Statistics &statistics)
{
    if (!good_weight_system(ws))
        return;

    sort(ws);
    weight_systems.insert(ws);
    ++statistics.weight_systems_found;
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
    Point x = history.points[n];

    for (int i = 0; i < n; ++i) {
        Long rel = history.point_weight_system_distances[i][i] -
                   history.point_weight_system_distances[n][i];
        if (rel > 0)
            return true;
        if (!disable_lex_compare && rel == 0 && history.points[i] < x)
            return true;
    }

    for (int i = 0; i < n; ++i) {
        Point x_other = history.points[i];
        Point x_diff = x_other - x;

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

    add_maybe(weight_systems, ws, statistics);

    switch (n) {
    case dim - 2:
        // The following happens when redundancies are not checked
        assert(builder.generator_count() == 2);
        // TODO: special case for last iteration brings a big performance boost
        if (defer_last_recursion) {
            pairs.insert(canonicalize(builder.to_pair()));
            ++statistics.final_pairs_found;

            // if (final_pairs_found % 10000 == 0)
            //     System.out.printf("%7.2f: builders: %d, unique: %d\n",
            //                       stopwatch.count(), buildersAddedCount,
            //                       builders.size());
            return;
        }
        break;
    case dim - 1:
        return;
    }

    auto symmetries = points_symmetries(span<Point>(history.points.data(), n));

    auto points = WeightSystemPointsBelow(ws);
    while (points.find_next()) {
        const Point &x = points.get();

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
        cout << stopwatch
             << " - weight systems: " << statistics.weight_systems_found
             << ", unique: " << weight_systems.size() << endl;
        cout << stopwatch
             << " - pairs: " << statistics.final_pairs_found
             << ", unique: " << pairs.size() << endl;

        for (const auto &pair : pairs) {
            WeightSystem ws = average(pair);

            // auto symmetries = points_symmetries(span<Point>(history.points.data(), n));

            auto points = WeightSystemPointsBelow(ws);
            while (points.find_next()) {
                const Point &x = points.get();

                if (!leads_to_allowed_weightsystem(x)
                    // || (!debug_ignore_symmetries && !is_sorted(x, symmetries))
                    )
                    continue;

                WeightSystem final_ws{};
                if (restrict(pair, x, final_ws))
                    add_maybe(weight_systems, final_ws, statistics);
            }
        }
    }

    cout << stopwatch
         << " - weight systems: " << statistics.weight_systems_found
         << ", unique: " << weight_systems.size() << endl;

    return 0;
}
