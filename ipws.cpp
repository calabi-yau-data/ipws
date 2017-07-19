#include <tclap/CmdLine.h>
#include <boost/filesystem.hpp>
#include <boost/optional.hpp>
#include <iostream>
#include <set>
#include <sstream>
#include <unordered_set>
#include "config.h"
#include "file.h"
#include "history.h"
#include "point.h"
#include "settings.h"
#include "stl_utils.h"
#include "stopwatch.h"
#include "weight_system.h"
#include "weight_system_builder.h"
#include "weight_system_pair.h"

using boost::optional;
using boost::none;
using gsl::span;
using std::array;
using std::cerr;
using std::cout;
using std::endl;
using std::unordered_set;
using std::string;
using std::vector;

namespace fs = boost::filesystem;

struct Statistics {
    unsigned weight_systems_found;
    unsigned final_pairs_found;
    unsigned ip_weight_systems;
};

void print_with_denominator(const WeightSystem &ws)
{
    Ring n = norm(ws);
    cout << n * r_denominator;
    for (const auto &w : ws.weights)
        cout << " " << w * r_numerator;
    cout << endl;
}

bool add_maybe(unordered_set<WeightSystem> &weight_systems, WeightSystem ws,
               Statistics &statistics)
{
    if (!good_weight_system(ws))
        return false;

    sort(ws);
    weight_systems.insert(ws);
    ++statistics.weight_systems_found;

    return true;
}

void rec(const WeightSystemBuilder &builder,
         unordered_set<WeightSystem> *weight_systems,
         unordered_set<WeightSystemPair> &final_pairs, unsigned n,
         History &history, Statistics &statistics, const Stopwatch &stopwatch)
{
    WeightSystem ws{};
    if (!builder.average_if_nonzero(ws))
        return;

    history.weight_systems[n] = ws;

    if (n < dim - g_settings.redundancy_check_skip_recursions &&
        last_point_redundant2(builder, n, history))
        return;

    if (weight_systems)
        add_maybe(*weight_systems, ws, statistics);

    switch (n) {
    case dim - 2:
        /* if (defer_last_recursion) */ {
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
            (!g_settings.debug_ignore_symmetries && !is_sorted(x, sym)))
            continue;

        history.points[n] = x;
        for (unsigned i = 0; i < n + 1; ++i)
            history.point_weight_system_distances[n][i] =
                distance(history.weight_systems[i], x);

        if (last_point_redundant(n, history))
            continue;

        rec(builder.restrict(x), weight_systems, final_pairs, n + 1, history,
            statistics, stopwatch);
    }
}

void write_sorted(File &f, const unordered_set<WeightSystem> &weight_systems)
{
    vector<WeightSystem> ws_list{};
    ws_list.reserve(weight_systems.size());
    std::copy(weight_systems.begin(), weight_systems.end(),
              std::back_inserter(ws_list));

    f.write(static_cast<uint32_t>(weight_systems.size()));
    for (const auto &ws : ws_list)
        write(f, ws);
}

void find_weight_systems_from_pair(unordered_set<WeightSystem> &weight_systems,
                                   const WeightSystemPair &pair,
                                   Statistics &statistics,
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
                if (g_settings.print_last_recursion_statistics) {
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

    if (g_settings.print_last_recursion_statistics) {
        size_t unique_count = weight_systems.size() - old_unique_count;
        cout << candidate_count << " " << current_weight_systems.size() << " "
             << unique_count << endl;
    }
}

void read_pairs(unordered_set<WeightSystem> *weight_systems,
                Statistics statistics, Stopwatch &stopwatch, File &pairs_in,
                const optional<fs::path> &write_weight_systems_from_pairs_dir)
{
    uint32_t size;
    pairs_in.read(size);

    srand(1234);
    for (unsigned i = 0; i < size; ++i) {
        // unsigned ws_size = 2 * dim * (sizeof uint16_t);
        // pairs_in.seek((rand() % (size / ws_size)) * ws_size);
        WeightSystemPair pair;

        read(pairs_in, pair.first);
        read(pairs_in, pair.second);

        unordered_set<WeightSystem> weight_systems_of_pair{};
        find_weight_systems_from_pair(weight_systems_of_pair, pair, statistics,
                                      stopwatch);

        if (write_weight_systems_from_pairs_dir) {
            // TODO
        }

        if (weight_systems)
            std::copy(weight_systems_of_pair.begin(),
                      weight_systems_of_pair.end(),
                      std::inserter(*weight_systems, weight_systems->end()));
    }
}

void classify(optional<File> &pairs_in, optional<File> &pairs_out,
              optional<File> &intermediate_ws_out,
              const optional<fs::path> &write_weight_systems_from_pairs_dir)
{
    Stopwatch stopwatch{};
    History history{};
    unordered_set<WeightSystem> weight_systems{};
    unordered_set<WeightSystemPair> final_pairs{};
    Statistics statistics{};

    if (pairs_in) {
        if (g_settings.ip_check || g_settings.count_weight_systems ||
            g_settings.print_candidates) {
            read_pairs(&weight_systems, statistics, stopwatch, *pairs_in,
                       write_weight_systems_from_pairs_dir);

            cerr << stopwatch
                 << " - weight systems: " << statistics.weight_systems_found
                 << ", unique: " << weight_systems.size() << endl;
        } else {
            read_pairs(nullptr, statistics, stopwatch, *pairs_in,
                       write_weight_systems_from_pairs_dir);
        }
    } else {
        if (intermediate_ws_out || g_settings.ip_check ||
            g_settings.count_weight_systems || g_settings.print_candidates) {
            rec(WeightSystemBuilder{}, &weight_systems, final_pairs, 0, history,
                statistics, stopwatch);

            cerr << stopwatch << " - intermediate weight systems: "
                 << statistics.weight_systems_found
                 << ", unique: " << weight_systems.size() << endl;
        } else {
            rec(WeightSystemBuilder{}, nullptr, final_pairs, 0, history,
                statistics, stopwatch);
        }

        cerr << stopwatch
             << " - weight system pairs: " << statistics.final_pairs_found
             << ", unique: " << final_pairs.size() << endl;

        if (intermediate_ws_out) {
            cerr << stopwatch << " - writing intermediate weight systems\n";

            write_sorted(*intermediate_ws_out, weight_systems);
            intermediate_ws_out = none;

            cerr << stopwatch << " - writing complete\n";
        }

        if (pairs_out) {
            cerr << stopwatch << " - writing pairs\n";

            pairs_out->write(static_cast<uint32_t>(final_pairs.size()));
            for (auto &pair : final_pairs) {
                write(*pairs_out, pair.first);
                write(*pairs_out, pair.second);
            }

            pairs_out = none;
            cerr << stopwatch << " - writing complete\n";
        }

        if (g_settings.ip_check || g_settings.count_weight_systems) {
            for (const auto &pair : final_pairs)
                find_weight_systems_from_pair(weight_systems, pair, statistics,
                                              stopwatch);

            cerr << stopwatch
                 << " - weight systems: " << statistics.weight_systems_found
                 << ", unique: " << weight_systems.size() << endl;
        }
    }

    if (g_settings.print_candidates)
        for (const auto &ws : weight_systems)
            print_with_denominator(ws);

    if (g_settings.ip_check) {
        for (const auto &ws : weight_systems)
            if (has_ip(ws))
                ++statistics.ip_weight_systems;

        cerr << stopwatch
             << " - weight systems: " << statistics.weight_systems_found
             << ", unique: " << weight_systems.size()
             << ", ip: " << statistics.ip_weight_systems << endl;
    }
}

int main(int argc, char *argv[])
{
    std::ostringstream description;
    description << "Classify weight systems with d=" << dim
                << ", r=" << r_numerator << "/" << r_denominator
                << " using ring type '" << typeid(Ring).name() << "'";

    TCLAP::CmdLine cmd(description.str(), ' ', GIT_REVISION);

    TCLAP::ValueArg<string> read_pairs_arg( //
        "", "read-pairs", "Read weight system pairs from file", false, "",
        "file", cmd);
    TCLAP::ValueArg<string> write_pairs_arg( //
        "", "write-pairs", "Write weight system pairs to file", false, "",
        "file", cmd);
    TCLAP::ValueArg<string> write_intermediate_weight_systems_arg(
        "", "write-intermediate", "Write intermediate weight systems to file",
        false, "", "file", cmd);
    TCLAP::ValueArg<string> write_weight_systems_from_pairs_arg(
        "", "write-weight-systems-from-pairs",
        "Write weight system obtained from pairs to directory", false, "",
        "directory", cmd);
    TCLAP::ValueArg<unsigned> skip_redundancy_check_arg(
        "", "skip-redundancy-check",
        "Skip redundancy check for the number of recursions (default 2)", false,
        2, "number", cmd);
    TCLAP::SwitchArg allow_weight_one_arg(
        "", "allow-weight-one", "Allow weight systems containing a weight of 1",
        cmd);
    TCLAP::SwitchArg allow_weight_one_half_arg(
        "", "allow-weight-one-half",
        "Allow weight systems containing a weight of 1/2", cmd);
    TCLAP::SwitchArg allow_weights_sum_one_arg(
        "", "allow-weights-sum-one",
        "Allow weight systems containing two weights with a sum of 1", cmd);
    TCLAP::SwitchArg print_stats_arg(
        "", "print-stats", "Print statistics during the last recursion", cmd);
    TCLAP::SwitchArg print_candidates_arg(
        "", "print-candidates", "Print weight system candidates", cmd);
    TCLAP::SwitchArg ignore_symmetries_arg(
        "", "ignore-symmetries", "Ignore symmetries (for debugging)", cmd);
    TCLAP::SwitchArg no_lex_order_arg(
        "", "no-lex-order", "Disable lexicographic order (for debugging)", cmd);
    TCLAP::SwitchArg ip_check_arg("", "ip-check", "Check IP property", cmd);
    TCLAP::SwitchArg count_weight_systems_arg("", "count-weight-systems",
                                              "Count weight systems", cmd);

    cmd.parse(argc, argv);

    g_settings.redundancy_check_skip_recursions =
        skip_redundancy_check_arg.getValue();
    g_settings.allow_weight_one = allow_weight_one_arg.getValue();
    g_settings.allow_weight_one_half = allow_weight_one_half_arg.getValue();
    g_settings.allow_weights_sum_one = allow_weights_sum_one_arg.getValue();
    g_settings.print_last_recursion_statistics = print_stats_arg.getValue();
    g_settings.print_candidates = print_candidates_arg.getValue();
    g_settings.debug_ignore_symmetries = ignore_symmetries_arg.getValue();
    g_settings.debug_disable_lex_order = no_lex_order_arg.getValue();
    g_settings.ip_check = ip_check_arg.getValue();
    g_settings.count_weight_systems = count_weight_systems_arg.getValue();

    optional<File> pairs_in{};
    if (read_pairs_arg.isSet()) {
        pairs_in = File::open(read_pairs_arg.getValue());
        if (!pairs_in) {
            cerr << "Could not open file '" << read_pairs_arg.getValue()
                 << "'\n";
            return EXIT_FAILURE;
        }
    }

    optional<File> pairs_out{};
    if (write_pairs_arg.isSet()) {
        fs::path path = write_pairs_arg.getValue();

        pairs_out = File::create_new(path);
        if (!pairs_out) {
            cerr << "Could not create new file " << path << endl;
            return EXIT_FAILURE;
        }
    }

    optional<File> intermediate_ws_out{};
    if (write_intermediate_weight_systems_arg.isSet()) {
        fs::path path = write_intermediate_weight_systems_arg.getValue();

        intermediate_ws_out = File::create_new(path);
        if (!intermediate_ws_out) {
            cerr << "Could not create new file " << path << endl;
            return EXIT_FAILURE;
        }
    }

    optional<fs::path> write_weight_systems_from_pairs_dir{};
    if (write_weight_systems_from_pairs_arg.isSet())
        write_weight_systems_from_pairs_dir =
            write_weight_systems_from_pairs_arg.getValue();

    try {
        if (!intermediate_ws_out && !g_settings.ip_check &&
            !g_settings.count_weight_systems && !pairs_out &&
            !write_weight_systems_from_pairs_dir &&
            !g_settings.print_candidates) {
            auto e = TCLAP::ArgException("No action specified");
            cmd.getOutput()->failure(cmd, e);
        }

        classify(pairs_in, pairs_out, intermediate_ws_out,
                 write_weight_systems_from_pairs_dir);
        return EXIT_SUCCESS;
    } catch (File::Error) {
        cerr << "IO error\n";
        return EXIT_FAILURE;
    } catch (TCLAP::ExitException) {
    }
}
