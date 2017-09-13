#include <tclap/CmdLine.h>
#include <algorithm>
#include <boost/optional.hpp>
#include <fstream>
#include <iostream>
#include <memory>
#include <random>
#include <set>
#include <sstream>
#include <unordered_set>
#include "buffered_reader.h"
#include "buffered_writer.h"
#include "config.h"
#include "history.h"
#include "point.h"
#include "settings.h"
#include "stl_utils.h"
#include "stopwatch.h"
#include "weight_system.h"
#include "weight_system_builder.h"
#include "weight_system_pair.h"

using boost::optional;
using gsl::span;
using std::cerr;
using std::cout;
using std::endl;
using std::unordered_set;
using std::string;
using std::vector;
using std::unique_ptr;
using std::make_unique;

struct Statistics {
    size_t weight_systems_found;
    size_t pairs_found;
    size_t ip_weight_systems;
};

void print_with_denominator(const WeightSystem<dim> &ws)
{
    Ring n = norm(ws);
    cout << n * r_denominator;
    for (const auto &w : ws.weights)
        cout << " " << w * r_numerator;
    cout << endl;
}

bool add_maybe(unordered_set<WeightSystem<dim>> &weight_systems,
               WeightSystem<dim> ws, Statistics &statistics)
{
    if (!good_weight_system(ws))
        return false;

    sort(ws);
    weight_systems.insert(ws);
    ++statistics.weight_systems_found;

    return true;
}

void rec(const WeightSystemBuilder &builder,
         unordered_set<WeightSystem<dim>> *weight_systems,
         unordered_set<WeightSystemPair> &pairs, unsigned n, History &history,
         Statistics &statistics, const Stopwatch &stopwatch)
{
    WeightSystem<dim> ws{};
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
        // defer last recursion
        pairs.insert(canonicalize(builder.to_pair()));
        ++statistics.pairs_found;

        if (statistics.pairs_found % 10000 == 0)
            cerr << stopwatch << " - pairs: " << statistics.pairs_found
                 << ", unique: " << pairs.size() << endl;

        return;
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

void write_config(BufferedWriter &f)
{
    write(f, static_cast<uint32_t>(dim));
    write(f, static_cast<uint32_t>(r_numerator));
    write(f, static_cast<uint32_t>(r_denominator));
}

void check_config(BufferedReader &f)
{
    uint32_t v;

    read(f, v);
    assert(v == dim);
    read(f, v);
    assert(v == r_numerator);
    read(f, v);
    assert(v == r_denominator);
}

void write_sorted(BufferedWriter &f,
                  const unordered_set<WeightSystem<dim>> &weight_systems)
{
    vector<WeightSystem<dim>> sorted{};
    sorted.reserve(weight_systems.size());

    std::copy(weight_systems.begin(), weight_systems.end(),
              std::back_inserter(sorted));

    std::sort(sorted.begin(), sorted.end());

    write_config(f);
    write(f, static_cast<uint64_t>(sorted.size()));
    for (const auto &ws : sorted)
        write_varint(f, ws);
}

void write_randomized(BufferedWriter &f,
                      const unordered_set<WeightSystemPair> &pairs)
{
    vector<WeightSystemPair> sorted{};
    sorted.reserve(pairs.size());

    std::copy(pairs.begin(), pairs.end(), std::back_inserter(sorted));

    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(sorted.begin(), sorted.end(), g);

    write_config(f);
    write(f, static_cast<uint64_t>(sorted.size()));
    for (const auto &pair : sorted) {
        write(f, pair.first);
        write(f, pair.second);
    }
}

void weight_systems_from_pair(unordered_set<WeightSystem<dim>> &weight_systems,
                              const WeightSystemPair &pair,
                              Statistics &statistics)
{
    WeightSystem<dim> ws = average(pair);

    // auto sym = symmetries(pair);

    auto points = WeightSystemPointsBelow(ws);
    while (points.find_next()) {
        const Point &x = points.get();

        if (!leads_to_allowed_weightsystem(x)
            // || (!debug_ignore_symmetries && !is_sorted(x, sym)))
            )
            continue;

        WeightSystem<dim> final_ws{};
        if (restrict(pair, x, final_ws))
            add_maybe(weight_systems, final_ws, statistics);
    }
}

void find_weight_systems(bool ip_only)
{
    Stopwatch stopwatch{};
    History history{};
    Statistics statistics{};
    unordered_set<WeightSystem<dim>> weight_systems{};
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

void find_pairs(BufferedWriter *ws_out, BufferedWriter *pairs_out)
{
    Stopwatch stopwatch{};
    History history{};
    Statistics statistics{};
    unordered_set<WeightSystem<dim>> weight_systems{};
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
        ws_out->flush();
        cerr << stopwatch << " - writing complete\n";
    }

    if (pairs_out) {
        cerr << stopwatch << " - writing pairs\n";
        write_randomized(*pairs_out, pairs);
        pairs_out->flush();
        cerr << stopwatch << " - writing complete\n";
    }
}

void find_weight_systems_from_pairs(BufferedReader &pairs_in, size_t start,
                                    optional<size_t> count_opt,
                                    BufferedWriter *ws_out)
{
    constexpr size_t pair_storage_size = 2 * WeightSystem<dim>::storage_size;
    unordered_set<WeightSystem<dim>> weight_systems{};
    Statistics statistics{};
    Stopwatch stopwatch{};
    unordered_set<WeightSystemPair> pairs{};
    // size_t candidate_count = 0;

    check_config(pairs_in);

    uint64_t pair_count;
    read(pairs_in, pair_count);

    size_t count = count_opt ? *count_opt : pair_count - start;

    cerr << stopwatch << " - total pairs: " << pair_count
         << ", pairs used: " << count << endl;

    assert(start + count <= pair_count);
    pairs_in.seek_relative(pair_storage_size * start);

    // srand(1234);
    for (size_t i = 0; i < count; ++i) {
        // pairs_in.seek_relative((rand() % (pair_count / pair_storage_size)) *
        //                        pair_storage_size);

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
        ws_out->flush();
        cerr << stopwatch << " - writing complete\n";
    }

    if (g_settings.print_weight_systems)
        for (const auto &ws : weight_systems)
            print_with_denominator(ws);
}

void check_ip(BufferedReader &ws_in)
{
    Stopwatch stopwatch{};

    check_config(ws_in);

    uint64_t ws_count;
    read(ws_in, ws_count);

    size_t ip_count = 0;
    for (size_t i = 0; i < ws_count; ++i) {
        WeightSystem<dim> ws{};
        read_varint(ws_in, ws);

        if (has_ip(ws)) {
            ++ip_count;
            if (g_settings.print_weight_systems)
                print_with_denominator(ws);
        }
    }

    cerr << stopwatch << " - weight systems: " << ws_count
         << ", ip: " << ip_count << endl;
}

void combine_ws_files(span<BufferedReader> ins, BufferedWriter &out)
{
    Stopwatch stopwatch{};

    vector<size_t> counts{};
    size_t total = 0;

    for (auto &in : ins) {
        check_config(in);

        uint64_t count;
        read(in, count);
        counts.push_back(count);
        total += count;

        cerr << "input weight systems: " << count << endl;
    }

    cerr << stopwatch << " - total input weight systems: " << total << endl;

    write_config(out);

    uint64_t count = 0;
    write(out, count);

    vector<WeightSystem<dim>> ws{};
    ws.resize(ins.size());

    for (size_t i = 0; i < counts.size(); ++i)
        if (counts[i] > 0)
            read_varint(ins[i], ws[i]);

    WeightSystem<dim> prev;

    while (true) {
        WeightSystem<dim> *smallest = nullptr;
        size_t smallest_index = 0;
        for (size_t i = 0; i < counts.size(); ++i) {
            if (counts[i] == 0)
                continue;

            if (smallest && *smallest < ws[i])
                continue;

            smallest = &ws[i];
            smallest_index = i;
        }
        if (!smallest)
            break;

        if (count == 0 || prev != *smallest) {
            ++count;
            write_varint(out, *smallest);
            prev = *smallest;
        }

        if (--counts[smallest_index] > 0)
            read_varint(ins[smallest_index], ws[smallest_index]);
    }

    // Write header again, now with correct count
    out.seek(0);
    write_config(out);
    write(out, count);

    cerr << stopwatch << " - combined weight systems: " << count << endl;
}

void diff_ws_files(BufferedReader &in1, BufferedReader &in2)
{
    Stopwatch stopwatch{};

    vector<size_t> counts{};

    check_config(in1);
    check_config(in2);

    uint64_t count1;
    read(in1, count1);
    uint64_t count2;
    read(in2, count2);

    cerr << stopwatch << " - input weight systems: " << count1 << " and "
         << count2 << endl;

    WeightSystem<dim> ws1{};
    WeightSystem<dim> ws2{};

    if (count1 > 0)
        read_varint(in1, ws1);
    if (count2 > 0)
        read_varint(in2, ws2);

    while (count1 > 0 && count2 > 0) {
        if (ws1 < ws2) {
            cout << "- " << ws1 << endl;
            if (--count1 > 0)
                read_varint(in1, ws1);
        } else if (ws2 < ws1) {
            cout << "+ " << ws2 << endl;
            if (--count2 > 0)
                read_varint(in2, ws2);
        } else {
            if (--count1 > 0)
                read_varint(in1, ws1);
            if (--count2 > 0)
                read_varint(in2, ws2);
        }
    }

    while (count1 > 0) {
        cout << "- " << ws1 << endl;
        if (--count1 > 0)
            read_varint(in1, ws1);
    }

    while (count2 > 0) {
        cout << "+ " << ws2 << endl;
        if (--count2 > 0)
            read_varint(in2, ws2);
    }
}

void print_count(BufferedReader &in)
{
    check_config(in);
    uint64_t count;
    read(in, count);
    cout << count << endl;
}

bool run(int argc, char *argv[])
{
    using std::fstream;
    using TCLAP::Arg;
    using TCLAP::SwitchArg;
    using TCLAP::ValueArg;

    std::ostringstream description;
    description << "Classify weight systems with d=" << dim
                << ", r=" << r_numerator << "/" << r_denominator
                << " using ring '" << ring_name << "'";

    TCLAP::CmdLine cmd(description.str(), ' ', GIT_REVISION);

    SwitchArg find_candidates_arg( //
        "", "find-candidates",
        "Find and optionally print IP weight system candidates");
    SwitchArg find_ip_arg( //
        "", "find-ip", "Find and optionally print IP weight systems");
    SwitchArg find_pairs_arg(
        "", "find-pairs",
        "Find weight system pairs in the penultimate recursion");
    SwitchArg print_count_arg( //
        "", "print-count", "Print the number entries in the given file");
    SwitchArg combine_ws_arg( //
        "", "combine-ws", "Combine given weight system file with a second one");
    SwitchArg diff_ws_arg( //
        "", "diff-ws", "Diff two given weight system files");

    vector<Arg *> arg_list;
    arg_list.push_back(&find_candidates_arg);
    arg_list.push_back(&find_ip_arg);
    arg_list.push_back(&find_pairs_arg);
    arg_list.push_back(&print_count_arg);
    arg_list.push_back(&combine_ws_arg);
    arg_list.push_back(&diff_ws_arg);
    cmd.xorAdd(arg_list);

    SwitchArg print_ws_arg( //
        "", "print-ws", "Enable printing of weight systems", cmd);
    ValueArg<string> ws_out_arg( //
        "", "ws-out", "Weight systems destination file", false, "", "file",
        cmd);
    ValueArg<string> ws_in_arg( //
        "", "ws-in", "Weight systems source file", false, "", "file", cmd);
    ValueArg<string> pairs_in_arg( //
        "", "pairs-in", "Pairs source file", false, "", "file", cmd);
    ValueArg<string> pairs_out_arg( //
        "", "pairs-out", "Pairs destination file", false, "", "file", cmd);
    ValueArg<size_t> start_arg( //
        "", "start", "", false, 0, "number", cmd);
    ValueArg<ptrdiff_t> count_arg( //
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

    TCLAP::UnlabeledMultiArg<string> files_arg( //
        "", "Additional files", false, "file", cmd);

    cmd.parse(argc, argv);

    g_settings.redundancy_check_skip_recursions =
        skip_redundancy_check_arg.getValue();
    g_settings.allow_weights_one = allow_weights_one_arg.getValue();
    g_settings.allow_weights_one_half = allow_weights_one_half_arg.getValue();
    g_settings.allow_weights_sum_one = allow_weights_sum_one_arg.getValue();
    g_settings.print_weight_systems = print_ws_arg.getValue();
    g_settings.debug_ignore_symmetries = ignore_symmetries_arg.getValue();
    g_settings.debug_disable_lex_order = no_lex_order_arg.getValue();

    size_t start = start_arg.getValue();
    optional<size_t> count{};
    if (count_arg.getValue() >= 0)
        count = count_arg.getValue();

    if (find_candidates_arg.getValue()) {
        if (pairs_in_arg.isSet()) {
            BufferedReader pairs_in(pairs_in_arg.getValue());

            unique_ptr<BufferedWriter> ws_out{};
            if (ws_out_arg.isSet()) {
                ws_out = make_unique<BufferedWriter>(ws_out_arg.getValue());
            }

            find_weight_systems_from_pairs(pairs_in, start, count,
                                           ws_out.get());
        } else {
            find_weight_systems(false);
        }
    } else if (find_ip_arg.getValue()) {
        if (ws_in_arg.isSet()) {
            BufferedReader ws_in(ws_in_arg.getValue());
            check_ip(ws_in);
        } else {
            find_weight_systems(true);
        }
    } else if (find_pairs_arg.getValue()) {
        unique_ptr<BufferedWriter> ws_out{};
        if (ws_out_arg.isSet()) {
            ws_out = make_unique<BufferedWriter>(ws_out_arg.getValue());
        }

        unique_ptr<BufferedWriter> pairs_out{};
        if (pairs_out_arg.isSet()) {
            pairs_out = make_unique<BufferedWriter>(pairs_out_arg.getValue());
        }

        find_pairs(ws_out.get(), pairs_out.get());
    } else if (print_count_arg.getValue()) {
        if (ws_in_arg.isSet()) {
            BufferedReader in(ws_in_arg.getValue());
            print_count(in);
        }
        if (pairs_in_arg.isSet()) {
            BufferedReader in(pairs_in_arg.getValue());
            print_count(in);
        }
    } else if (combine_ws_arg.isSet()) {
        if (ws_out_arg.isSet()) {
            BufferedWriter ws_out(ws_out_arg.getValue());

            vector<BufferedReader> ws_in{};
            for (const auto &filename : files_arg.getValue())
                ws_in.push_back(BufferedReader{filename});

            combine_ws_files(span<BufferedReader>(ws_in), ws_out);
        }
    } else if (diff_ws_arg.isSet() && files_arg.getValue().size() == 2) {
        BufferedReader ws1_in{files_arg.getValue()[0]};
        BufferedReader ws2_in{files_arg.getValue()[1]};

        diff_ws_files(ws1_in, ws2_in);
    }
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
