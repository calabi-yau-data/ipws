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
    unsigned pairs_found;
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
         unordered_set<WeightSystemPair> &pairs, unsigned n, History &history,
         Statistics &statistics, const Stopwatch &stopwatch)
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
            pairs.insert(canonicalize(builder.to_pair()));
            ++statistics.pairs_found;

            if (statistics.pairs_found % 10000 == 0)
                cerr << stopwatch << " - pairs: " << statistics.pairs_found
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
            (!g_settings.debug_ignore_symmetries && !is_sorted(x, sym)))
            continue;

        history.points[n] = x;
        for (unsigned i = 0; i < n + 1; ++i)
            history.point_weight_system_distances[n][i] =
                distance(history.weight_systems[i], x);

        if (last_point_redundant(n, history))
            continue;

        rec(builder.restrict(x), weight_systems, pairs, n + 1, history,
            statistics, stopwatch);
    }
}

void write_config(File &f)
{
    f.write(static_cast<uint32_t>(dim));
    f.write(static_cast<uint32_t>(r_numerator));
    f.write(static_cast<uint32_t>(r_denominator));
}

void check_config(File &f)
{
    uint32_t v;

    f.read(v);
    assert(v == dim);
    f.read(v);
    assert(v == r_numerator);
    f.read(v);
    assert(v == r_denominator);
}

void write_sorted(File &f, const unordered_set<WeightSystem> &weight_systems)
{
    vector<WeightSystem> ws_list{};
    ws_list.reserve(weight_systems.size());
    std::copy(weight_systems.begin(), weight_systems.end(),
              std::back_inserter(ws_list));

    write_config(f);

    f.write(static_cast<uint32_t>(weight_systems.size()));
    for (const auto &ws : ws_list)
        write(f, ws);
}

void weight_systems_from_pair(unordered_set<WeightSystem> &weight_systems,
                              const WeightSystemPair &pair,
                              Statistics &statistics)
{
    WeightSystem ws = average(pair);

    // auto sym = symmetries(pair);

    auto points = WeightSystemPointsBelow(ws);
    while (points.find_next()) {
        const Point &x = points.get();

        if (!leads_to_allowed_weightsystem(x)
            // || (!debug_ignore_symmetries && !is_sorted(x, sym)))
            )
            continue;

        WeightSystem final_ws{};
        if (restrict(pair, x, final_ws))
            add_maybe(weight_systems, final_ws, statistics);
    }
}

void find_weight_systems(bool ip_only)
{
    Stopwatch stopwatch{};
    History history{};
    Statistics statistics{};
    unordered_set<WeightSystem> weight_systems{};
    unordered_set<WeightSystemPair> pairs{};

    rec(WeightSystemBuilder{}, &weight_systems, pairs, 0, history, statistics,
        stopwatch);
    cerr << stopwatch << " - intermediate weight systems: "
         << statistics.weight_systems_found
         << ", unique: " << weight_systems.size() << endl;

    for (const auto &pair : pairs)
        weight_systems_from_pair(weight_systems, pair, statistics);

    cerr << stopwatch
         << " - weight systems: " << statistics.weight_systems_found
         << ", unique: " << weight_systems.size() << endl;

    if (ip_only) {
        for (const auto &ws : weight_systems) {
            if (has_ip(ws)) {
                ++statistics.ip_weight_systems;
                if (g_settings.print_weight_systems)
                    print_with_denominator(ws);
            }
        }
        cerr << stopwatch
             << " - weight systems: " << statistics.weight_systems_found
             << ", unique: " << weight_systems.size()
             << ", ip: " << statistics.ip_weight_systems << endl;
    } else if (g_settings.print_weight_systems) {
        for (const auto &ws : weight_systems)
            print_with_denominator(ws);
    }
}

void find_pairs(optional<File> &ws_out, optional<File> &pairs_out)
{
    Stopwatch stopwatch{};
    History history{};
    Statistics statistics{};
    unordered_set<WeightSystem> weight_systems{};
    unordered_set<WeightSystemPair> pairs{};

    rec(WeightSystemBuilder{}, &weight_systems, pairs, 0, history, statistics,
        stopwatch);
    cerr << stopwatch << " - pairs: " << statistics.pairs_found
         << ", unique: " << pairs.size() << endl;
    cerr << stopwatch
         << " - weight systems: " << statistics.weight_systems_found
         << ", unique: " << weight_systems.size() << endl;

    if (ws_out) {
        cerr << stopwatch << " - writing weight systems\n";
        write_sorted(*ws_out, weight_systems);
        ws_out = none;
        cerr << stopwatch << " - writing complete\n";
    }

    if (pairs_out) {
        cerr << stopwatch << " - writing pairs\n";

        write_config(*pairs_out);

        pairs_out->write(static_cast<uint32_t>(pairs.size()));
        for (auto &pair : pairs) {
            write(*pairs_out, pair.first);
            write(*pairs_out, pair.second);
        }
        pairs_out = none;
        cerr << stopwatch << " - writing complete\n";
    }
}

void find_weight_systems_from_pairs(File &pairs_in, unsigned start,
                                    optional<unsigned> count_opt,
                                    optional<File> &ws_out)
{
    // TODO: Only generate weight systems that were not already written to file.
    // TODO: Do not create file in final location until complete.
    unordered_set<WeightSystem> weight_systems{};
    Statistics statistics{};
    Stopwatch stopwatch{};
    unordered_set<WeightSystemPair> pairs{};
    // size_t candidate_count = 0;

    check_config(pairs_in);

    uint32_t pair_count;
    pairs_in.read(pair_count);

    unsigned count = count_opt ? *count_opt : pair_count - start;

    cerr << stopwatch << " - total pairs: " << pair_count
         << ", pairs used: " << pair_count << endl;

    assert(start + count <= pair_count);
    pairs_in.seek_relative(weight_system_storage_size * start);

    for (unsigned i = 0; i < count; ++i) {
        WeightSystemPair pair{};

        read(pairs_in, pair.first);
        read(pairs_in, pair.second);

        weight_systems_from_pair(weight_systems, pair, statistics);

        // size_t new_candidate_count = candidate_count + weight_systems.size();
        // if (candidate_count / 100000 != new_candidate_count / 100000)
        //     cerr << stopwatch << " - weight systems: " << new_candidate_count
        //          << endl;
        // candidate_count = new_candidate_count;
    }

    cerr << stopwatch
         << " - weight systems: " << statistics.weight_systems_found
         << ", unique: " << weight_systems.size() << endl;

    if (ws_out) {
        cerr << stopwatch << " - writing weight systems\n";
        write_sorted(*ws_out, weight_systems);
        ws_out = none;
        cerr << stopwatch << " - writing complete\n";
    }

    if (g_settings.print_weight_systems)
        for (const auto &ws : weight_systems)
            print_with_denominator(ws);
}

void check_ip(File &ws_in)
{
    Stopwatch stopwatch{};

    check_config(ws_in);

    uint32_t ws_count;
    ws_in.read(ws_count);

    unsigned ip_count = 0;
    for (unsigned i = 0; i < ws_count; ++i) {
        WeightSystem ws{};
        read(ws_in, ws);

        if (has_ip(ws)) {
            ++ip_count;
            if (g_settings.print_weight_systems)
                print_with_denominator(ws);
        }
    }

    cerr << stopwatch << " - weight systems: " << ws_count
         << ", ip: " << ip_count << endl;
}

void combine_ws_files(File &in1, File &in2, File &out)
{
    check_config(in1);
    check_config(in2);

    uint32_t count1;
    uint32_t count2;
    in1.read(count1);
    in2.read(count2);

    write_config(out);

    WeightSystem ws1{};
    WeightSystem ws2{};

    if (count1 > 0)
        read(in1, ws1);
    if (count2 > 0)
        read(in2, ws2);

    while (count1 > 0 && count2 > 0) {
        if (ws1 < ws2) {
            write(out, ws1);
            if (--count1 > 0)
                read(in1, ws1);
        } else if (ws2 < ws1) {
            write(out, ws2);
            if (--count2 > 0)
                read(in2, ws2);
        } else {
            if (--count1 > 0)
                read(in1, ws1);
        }
    }

    while (count1 > 0) {
        write(out, ws1);
        if (--count1 > 0)
            read(in1, ws1);
    }

    while (count2 > 0) {
        write(out, ws2);
        if (--count2 > 0)
            read(in2, ws2);
    }
}

File open_file(const fs::path &path)
{
    optional<File> f = File::open(path);
    if (!f) {
        cerr << "Could not open file " << path << endl;
        exit(EXIT_FAILURE);
    }
    return *f;
}

File create_file(const fs::path &path)
{
    optional<File> f = File::create_new(path);
    if (!f) {
        cerr << "Could not create new file " << path << endl;
        exit(EXIT_FAILURE);
    }
    return *f;
}

bool run(int argc, char *argv[])
{
    using TCLAP::Arg;
    using TCLAP::SwitchArg;
    using TCLAP::ValueArg;

    std::ostringstream description;
    description << "Classify weight systems with d=" << dim
                << ", r=" << r_numerator << "/" << r_denominator
                << " using ring type '" << typeid(Ring).name() << "'";

    TCLAP::CmdLine cmd(description.str(), ' ', GIT_REVISION);

    SwitchArg find_candidates_arg( //
        "", "find-candidates",
        "Find and optionally print IP weight system candidates");
    SwitchArg find_ip_arg( //
        "", "find-ip", "Find and optionally print IP weight systems");
    SwitchArg find_pairs_arg(
        "", "find-pairs",
        "Find weight system pairs in the penultimate recursion");
    // SwitchArg combine_ws_arg("", "combine-ws",
    //                          "Combine to weight systems files");

    vector<Arg *> arg_list;
    arg_list.push_back(&find_candidates_arg);
    arg_list.push_back(&find_ip_arg);
    arg_list.push_back(&find_pairs_arg);
    // arg_list.push_back(&combine_ws_arg);
    cmd.xorAdd(arg_list);

    SwitchArg print_ws_arg( //
        "", "print-ws", "Enable printing of weight systems", cmd);
    ValueArg<string> ws_out_arg( //
        "", "ws-out", "Weight systems destination file", false, "", "file",
        cmd);
    ValueArg<string> ws_in_arg( //
        "", "ws-in", "Weight systems source file", false, "", "file", cmd);
    ValueArg<string> ws_in2_arg( //
        "", "ws-in2", "Weight systems source file", false, "", "file", cmd);
    ValueArg<string> pairs_in_arg( //
        "", "pairs-in", "Pairs source file", false, "", "file", cmd);
    ValueArg<string> pairs_out_arg( //
        "", "pairs-out", "Pairs destination file", false, "", "file", cmd);
    ValueArg<unsigned> start_arg( //
        "", "start", "", false, 0, "number", cmd);
    ValueArg<int> count_arg( //
        "", "count", "", false, -1, "number", cmd);

    ValueArg<unsigned> skip_redundancy_check_arg(
        "", "skip-redundancy-check",
        "Skip redundancy check for the number of recursions (default 2)", false,
        2, "number", cmd);
    SwitchArg allow_weights_one_arg(
        "", "allow-weights-one",
        "Allow weight systems containing a weight of 1", cmd);
    SwitchArg allow_weights_one_half_arg(
        "", "allow-weights-one-half",
        "Allow weight systems containing a weight of 1/2", cmd);
    SwitchArg allow_weights_sum_one_arg(
        "", "allow-weights-sum-one",
        "Allow weight systems containing two weights with a sum of 1", cmd);
    SwitchArg ignore_symmetries_arg("", "ignore-symmetries",
                                    "Ignore symmetries (for debugging)", cmd);
    SwitchArg no_lex_order_arg(
        "", "no-lex-order", "Disable lexicographic order (for debugging)", cmd);

    cmd.parse(argc, argv);

    g_settings.redundancy_check_skip_recursions =
        skip_redundancy_check_arg.getValue();
    g_settings.allow_weights_one = allow_weights_one_arg.getValue();
    g_settings.allow_weights_one_half = allow_weights_one_half_arg.getValue();
    g_settings.allow_weights_sum_one = allow_weights_sum_one_arg.getValue();
    g_settings.print_weight_systems = print_ws_arg.getValue();
    g_settings.debug_ignore_symmetries = ignore_symmetries_arg.getValue();
    g_settings.debug_disable_lex_order = no_lex_order_arg.getValue();

    unsigned start = start_arg.getValue();
    optional<unsigned> count{};
    if (count_arg.getValue() >= 0)
        count = count_arg.getValue();

    if (find_candidates_arg.getValue()) {
        if (pairs_in_arg.isSet()) {
            File pairs_in = open_file(pairs_in_arg.getValue());

            optional<File> ws_out{};
            if (ws_out_arg.isSet())
                ws_out = create_file(ws_out_arg.getValue());

            find_weight_systems_from_pairs(pairs_in, start, count, ws_out);
        } else {
            find_weight_systems(false);
        }
    } else if (find_ip_arg.getValue()) {
        if (ws_in_arg.isSet()) {
            File ws_in = open_file(ws_in_arg.getValue());
            check_ip(ws_in);
        } else {
            find_weight_systems(true);
        }
    } else if (find_pairs_arg.getValue()) {
        optional<File> ws_out{};
        if (ws_out_arg.isSet())
            ws_out = create_file(ws_out_arg.getValue());

        optional<File> pairs_out{};
        if (pairs_out_arg.isSet())
            pairs_out = create_file(pairs_out_arg.getValue());

        find_pairs(ws_out, pairs_out);
    }
    // else if (combine_ws_arg.getValue()) {
    //     combine_ws_files(*ws_in, *ws_in2, *ws_out);
    // }
    return true;
}

int main(int argc, char *argv[])
{
    try {
        return run(argc, argv) ? EXIT_SUCCESS : EXIT_FAILURE;
    } catch (const std::exception &e) {
        cerr << e.what() << endl;
        return EXIT_FAILURE;
    }
}
