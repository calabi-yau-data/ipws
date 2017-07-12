#include <arpa/inet.h>
#include <fstream>
#include <iostream>
#include <set>
#include <unordered_set>
#include "config.h"
#include "point.h"
#include "stl_utils.h"
#include "stopwatch.h"
#include "weight_system.h"
#include "weight_system_builder.h"
#include "weight_system_pair.h"

using gsl::span;
using std::array;
using std::cerr;
using std::cout;
using std::endl;
using std::pair;
using std::unordered_set;
using std::set;
using std::setprecision;
using std::vector;

set<WeightSystemPair> pairs; // TODO: unordered_set? TODO: not global

struct History {
    array<Point, dim - 1> points;
    array<WeightSystem, dim> weight_systems;
    array<array<Ring, dim - 1>, dim - 1> point_weight_system_distances;
};

struct Statistics {
    unsigned weight_systems_found;
    unsigned final_pairs_found;
    unsigned ip_weight_systems;
};

using WeightSystemCollection = unordered_set<WeightSystem>;

void print_with_denominator(const WeightSystem &ws)
{
    Ring n = norm(ws);
    cout << n * r_denominator;
    for (const auto &w : ws.weights)
        cout << " " << w * r_numerator;
    cout << endl;
}

bool good_weight_system(const WeightSystem &ws)
{
    Ring n = norm(ws);

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

bool add_maybe(WeightSystemCollection &weight_systems, WeightSystem ws,
               Statistics &statistics)
{
    if (!good_weight_system(ws))
        return false;

    sort(ws);
    weight_systems.insert(ws);
    ++statistics.weight_systems_found;

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
            Ring diff = history.point_weight_system_distances[i][i] -
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
        Ring rel = history.point_weight_system_distances[i][i] -
                   history.point_weight_system_distances[n][i];
        if (rel > 0)
            return true;
        if (!disable_lex_compare && rel == 0 && history.points[i] < x)
            return true;
    }

    for (int i = 0; i < n; ++i) {
        Point x_other = history.points[i];
        Point x_diff = x_other - x;

        Ring v = gcd(x_diff);
        if (v != 1)
            x_diff /= v;

        bool all_positive = true;
        for (unsigned j = 0; j < dim; ++j) {
            if (x_other.coords[j] + x_diff.coords[j] < 0) {
                all_positive = false;
                break;
            }
        }
        if (all_positive)
            return true;

        all_positive = true;
        for (unsigned j = 0; j < dim; ++j) {
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
         const WeightSystemBuilder &builder, unsigned n, History &history,
         Statistics &statistics, const Stopwatch &stopwatch)
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
        if (defer_last_recursion || write_weight_system_pairs) {
            pairs.insert(canonicalize(builder.to_pair()));
            ++statistics.final_pairs_found;

            if (statistics.final_pairs_found % 10000 == 0)
                cerr << stopwatch
                     << " - pairs: " << statistics.final_pairs_found
                     << ", unique: " << pairs.size() << endl;

            return;
        }
        break;
    case dim - 1:
        return;
    }

    auto sym = symmetries(span<const Point>(history.points.data(), n));

    auto points = WeightSystemPointsBelow(ws);
    while (points.find_next()) {
        const Point &x = points.get();

        if (!leads_to_allowed_weightsystem(x) ||
            (!debug_ignore_symmetries && !is_sorted(x, sym)))
            continue;

        history.points[n] = x;
        for (unsigned i = 0; i < n + 1; ++i)
            history.point_weight_system_distances[n][i] =
                distance(history.weight_systems[i], x);

        if (last_point_redundant(n, history))
            continue;

        rec(weight_systems, builder.restrict(x), n + 1, history, statistics,
            stopwatch);
    }
}

void process_pair(WeightSystemCollection &weight_systems,
                  const WeightSystemPair &pair, Statistics &statistics,
                  const Stopwatch &stopwatch)
{
    WeightSystem ws = average(pair);
    unsigned candidate_count = 0;
    size_t old_unique_count = weight_systems.size();
    unordered_set<WeightSystem> current_weight_systems{};

    // auto sym = symmetries(pair);

    auto points = WeightSystemPointsBelow(ws);
    while (points.find_next()) {
        const Point &x = points.get();

        if (!leads_to_allowed_weightsystem(x)
            // || (!debug_ignore_symmetries && !is_sorted(x, sym)))
            )
            continue;

        WeightSystem final_ws{};
        if (restrict(pair, x, final_ws)) {
            if (add_maybe(weight_systems, final_ws, statistics)) {
                if (print_last_recursion_statistics) {
                    ++candidate_count;
                    current_weight_systems.insert(final_ws);
                }

                if (statistics.weight_systems_found % 1000000 == 0)
                    cerr << stopwatch << " - weight systems: "
                         << statistics.weight_systems_found
                         << ", unique: " << weight_systems.size() << endl;
            }
        }
    }

    if (print_last_recursion_statistics) {
        size_t unique_count = weight_systems.size() - old_unique_count;
        cout << candidate_count << " " << current_weight_systems.size() << " "
             << unique_count << endl;
    }
}

int main()
{
    Stopwatch stopwatch{};

    History history{};
    WeightSystemCollection weight_systems{};
    Statistics statistics{};

    if (read_weight_system_pairs) {
        std::ifstream cones_in{"pairs", std::ifstream::binary};

        srand(1234);
        while (cones_in) {
            // cones_in.seekg((rand() % 46739902) * 20);
            WeightSystemPair pair;

            for (unsigned i = 0; i < dim; ++i) {
                uint16_t v16;
                cones_in.read(reinterpret_cast<char *>(&v16), sizeof(v16));
                pair.first.weights[i] = ntohs(v16);
            }

            for (unsigned i = 0; i < dim; ++i) {
                uint16_t v16;
                cones_in.read(reinterpret_cast<char *>(&v16), sizeof(v16));
                pair.second.weights[i] = ntohs(v16);
            }

            if (!cones_in)
                break;

            process_pair(weight_systems, pair, statistics, stopwatch);
        }
    } else {
        rec(weight_systems, WeightSystemBuilder{}, 0, history, statistics,
            stopwatch);
    }

    if (write_weight_system_pairs) {
        cerr << stopwatch << " - writing\n";

        std::ofstream pairs_out{"pairs", std::ofstream::binary};
        for (auto &pair : pairs) {
            for (unsigned i = 0; i < dim; ++i) {
                auto v = pair.first.weights[i];
                assert(v >= 0 && v <= UINT16_MAX);
                uint16_t v16 = htons(static_cast<uint16_t>(v));
                pairs_out.write(reinterpret_cast<const char *>(&v16),
                                sizeof(v16));
            }
            for (unsigned i = 0; i < dim; ++i) {
                auto v = pair.second.weights[i];
                assert(v >= 0 && v <= UINT16_MAX);
                uint16_t v16 = htons(static_cast<uint16_t>(v));
                pairs_out.write(reinterpret_cast<const char *>(&v16),
                                sizeof(v16));
            }
        }
    }

    if (defer_last_recursion) {
        cerr << stopwatch
             << " - weight systems: " << statistics.weight_systems_found
             << ", unique: " << weight_systems.size() << endl;
        cerr << stopwatch << " - pairs: " << statistics.final_pairs_found
             << ", unique: " << pairs.size() << endl;

        for (const auto &pair : pairs)
            process_pair(weight_systems, pair, statistics, stopwatch);
    }

    cerr << stopwatch
         << " - weight systems: " << statistics.weight_systems_found
         << ", unique: " << weight_systems.size() << endl;

    for (const auto &ws : weight_systems) {
        if (print_candidates)
            print_with_denominator(ws);
        if (has_ip(ws))
            ++statistics.ip_weight_systems;
    }

    cerr << stopwatch
         << " - weight systems: " << statistics.weight_systems_found
         << ", unique: " << weight_systems.size()
         << ", ip: " << statistics.ip_weight_systems << endl;

    return 0;
}
