#include <argagg/argagg.hpp>
#include <experimental/optional>
#include <iostream>
#include <set>
#include <unordered_set>
#include "config.h"
#include "file.h"
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
using std::experimental::optional;
using std::pair;
using std::unordered_set;
using std::set;
using std::setprecision;
using std::string;
using std::vector;

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

struct Settings {
    unsigned redundancy_check_skip_recursions;
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
            if (diff > 0 || (!debug_disable_lex_compare && diff == 0 &&
                             x > history.points[i]))
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
        if (!debug_disable_lex_compare && rel == 0 && history.points[i] < x)
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

void rec(const WeightSystemBuilder &builder,
         WeightSystemCollection &weight_systems,
         set<WeightSystemPair> &final_pairs, unsigned n, History &history,
         Statistics &statistics, const Stopwatch &stopwatch,
         const Settings &settings)
{
    WeightSystem ws{};
    if (!builder.average_if_nonzero(ws))
        return;

    history.weight_systems[n] = ws;

    if (n < dim - settings.redundancy_check_skip_recursions &&
        last_point_redundant2(builder, n, history))
        return;

    add_maybe(weight_systems, ws, statistics);

    switch (n) {
    case dim - 2:
        // The following happens when redundancies are not checked
        assert(builder.generator_count() == 2);

        if (defer_last_recursion) {
            final_pairs.insert(canonicalize(builder.to_pair()));
            ++statistics.final_pairs_found;

            if (statistics.final_pairs_found % 10000 == 0)
                cerr << stopwatch << " - weight system pairs: "
                     << statistics.final_pairs_found
                     << ", unique: " << final_pairs.size() << endl;

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

        rec(builder.restrict(x), weight_systems, final_pairs, n + 1, history,
            statistics, stopwatch, settings);
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

bool classify(optional<File> &pairs_in, optional<File> &pairs_out,
              const Settings &settings)
{
    Stopwatch stopwatch{};
    History history{};
    WeightSystemCollection weight_systems{};
    set<WeightSystemPair> final_pairs{}; // TODO: unordered_set?
    Statistics statistics{};

    if (pairs_in) {
        try {
            uint32_t size;
            pairs_in->read(size);

            srand(1234);
            for (unsigned i = 0; i < size; ++i) {
                // unsigned ws_size = 2 * dim * (sizeof uint16_t);
                // pairs_in.seek((rand() % (size / ws_size)) * ws_size);
                WeightSystemPair pair;
                array<uint16_t, dim> data;

                pairs_in->read(span<uint16_t>(data));
                for (unsigned i = 0; i < dim; ++i)
                    pair.first.weights[i] = data[i];

                pairs_in->read(span<uint16_t>(data));
                for (unsigned i = 0; i < dim; ++i)
                    pair.second.weights[i] = data[i];

                process_pair(weight_systems, pair, statistics, stopwatch);
            }
        } catch (File::Error) {
            cerr << "Read error\n";
            return false;
        }
    } else {
        rec(WeightSystemBuilder{}, weight_systems, final_pairs, 0, history,
            statistics, stopwatch, settings);
    }

    if (pairs_out) {
        try {
            cerr << stopwatch << " - writing\n";

            pairs_out->write(static_cast<uint32_t>(final_pairs.size()));

            for (auto &pair : final_pairs) {
                array<uint16_t, dim> data;

                for (unsigned i = 0; i < dim; ++i) {
                    auto v = pair.first.weights[i];
                    assert(v >= 0 && v <= UINT16_MAX);
                    data[i] = static_cast<uint16_t>(v);
                }
                pairs_out->write(span<uint16_t>(data));

                for (unsigned i = 0; i < dim; ++i) {
                    auto v = pair.second.weights[i];
                    assert(v >= 0 && v <= UINT16_MAX);
                    data[i] = static_cast<uint16_t>(v);
                }
                pairs_out->write(span<uint16_t>(data));
            }
        } catch (File::Error) {
            cerr << "Write error\n";
            return false;
        }
    }

    if (defer_last_recursion) {
        cerr << stopwatch
             << " - weight systems: " << statistics.weight_systems_found
             << ", unique: " << weight_systems.size() << endl;
        cerr << stopwatch
             << " - weight system pairs: " << statistics.final_pairs_found
             << ", unique: " << final_pairs.size() << endl;

        for (const auto &pair : final_pairs)
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

    return true;
}

int main(int argc, char *argv[])
{
    argagg::parser argparser{{
        {"help", {"-h", "--help"}, "Show this help message", 0},
        {
            "read-pairs",
            {"--read-pairs"},
            "Read weight system pairs from given file",
            1,
        },
        {
            "write-pairs",
            {"--write-pairs"},
            "Write weight system pairs to given file",
            1,
        },
        {
            "skip-redundancy-check",
            {"--skip-redundancy-check"},
            "Skip redundancy check on the number of recursions given "
            "(default 2)",
            1,
        },
    }};

    Settings settings{};
    settings.redundancy_check_skip_recursions = 2;

    argagg::parser_results args{};
    try {
        args = argparser.parse(argc, argv);

        if (args["skip-redundancy-check"])
            settings.redundancy_check_skip_recursions =
                args["skip-redundancy-check"];

    } catch (const std::exception &e) {
        cerr << e.what() << endl;
        return EXIT_FAILURE;
    }

    if (args["help"]) {
        cerr << "Classify weight systems with d=" << dim
             << ", r=" << r_numerator << "/" << r_denominator << endl
             << "Revision: " << GIT_REVISION << endl;

        if (defer_last_recursion)
            cerr << "Last recursion is deferred\n";

        if (allow_weight_one || allow_weight_one_half || allow_weights_sum_one)
            cerr << "### Additional weight systems are classified ###\n";

        if (debug_ignore_symmetries || debug_disable_lex_compare)
            cerr << "### Some features are disabled ### \n";

        cerr << argparser;
        return EXIT_SUCCESS;
    }

    auto read_pairs_path = args["read-pairs"].as<string>("");
    auto write_pairs_path = args["write-pairs"].as<string>("");

    optional<File> pairs_in{};
    if (!read_pairs_path.empty()) {
        pairs_in = File::open(read_pairs_path);
        if (!pairs_in) {
            cerr << "Could not open file '" << read_pairs_path << "'\n";
            return EXIT_FAILURE;
        }
    }

    optional<File> pairs_out{};
    if (!write_pairs_path.empty()) {
        if (!defer_last_recursion) {
            cerr << "Compile with defer_last_recursion = true!\n";
            return EXIT_FAILURE;
        }

        pairs_out = File::create_new(write_pairs_path);
        if (!pairs_out) {
            cerr << "Could not create new file '" << write_pairs_path << "'\n";
            return EXIT_FAILURE;
        }
    }

    return classify(pairs_in, pairs_out, settings) ? EXIT_SUCCESS
                                                   : EXIT_FAILURE;
}
