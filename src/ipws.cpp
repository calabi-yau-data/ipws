#include "buffered_reader.h"
#include "buffered_writer.h"
#include "config.h"
#include "history.h"
#include "point.h"
#include "polytope.h"
#include "read_write.h"
#include "settings.h"
#include "stl_utils.h"
#include "stopwatch.h"
#include "weight_system.h"
#include "weight_system_builder.h"
#include "weight_system_pair.h"
#include <algorithm>
#include <boost/optional.hpp>
#include <fstream>
#include <iostream>
#include <memory>
#include <random>
#include <set>
#include <sstream>
#include <tclap/CmdLine.h>
#include <unordered_set>

using boost::optional;
using gsl::span;
using std::array;
using std::cerr;
using std::cout;
using std::endl;
using std::make_unique;
using std::string;
using std::unique_ptr;
using std::unordered_set;
using std::vector;

struct Statistics {
    unsigned long weight_systems_found;
    unsigned long pairs_found;
    unsigned long ip_weight_systems;
};

void print_with_denominator(std::ostream &os, const WeightSystem<dim> &ws)
{
    Ring n = norm(ws);
    os << n * r_denominator;
    for (const auto &w : ws.weights)
        os << " " << w * r_numerator;
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

    // cout << builder;
    // cout << "q" << n << " = " << ws << endl;

    if (n + g_settings.redundancy_check_skip_recursions < dim &&
        last_point_redundant(builder, n, history))
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

        // cout << "x" << n + 1 << " = " << x << endl;

        history.points[n] = x;
        for (unsigned i = 0; i < n + 1; ++i)
            history.point_weight_system_distances[n][i] =
                distance(history.weight_systems[i], x);

        if (last_point_redundant_quick(n, history))
            continue;

        rec(builder.restrict(x), weight_systems, pairs, n + 1, history,
            statistics, stopwatch);
    }
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
    write64u(f, sorted.size());
    for (const auto &ws : sorted)
        write_varint(f, ws);
}

void write_randomized(BufferedWriter &f,
                      const unordered_set<WeightSystemPair> &pairs)
{
    vector<WeightSystemPair> sorted{};
    sorted.reserve(pairs.size());

    std::copy(pairs.begin(), pairs.end(), std::back_inserter(sorted));

    std::mt19937 g{};
    std::shuffle(sorted.begin(), sorted.end(), g);

    write_config(f);
    write64u(f, sorted.size());
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
                if (g_settings.print_weight_systems) {
                    print_with_denominator(cout, ws);
                    cout << endl;
                }
            }
        }
        cerr << stopwatch
             << " - weight systems: " << statistics.weight_systems_found
             << ", unique: " << weight_systems.size()
             << ", ip: " << statistics.ip_weight_systems << endl;
    } else if (g_settings.print_weight_systems) {
        for (const auto &ws : weight_systems) {
            print_with_denominator(cout, ws);
            cout << endl;
        }
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

void find_weight_systems_from_pairs(BufferedReader &pairs_in,
                                    unsigned long start,
                                    optional<unsigned long> count_opt,
                                    BufferedWriter *ws_out)
{
    constexpr size_t pair_storage_size = 2 * WeightSystem<dim>::storage_size;
    unordered_set<WeightSystem<dim>> weight_systems{};
    Statistics statistics{};
    Stopwatch stopwatch{};
    unordered_set<WeightSystemPair> pairs{};
    // unsigned long candidate_count = 0;

    check_config(pairs_in);

    uint64_t pair_count;
    read(pairs_in, pair_count);

    unsigned long count = count_opt ? *count_opt : pair_count - start;

    cerr << stopwatch << " - total pairs: " << pair_count
         << ", pairs used: " << count << endl;

    assert(start + count <= pair_count);
    pairs_in.seek_relative(pair_storage_size * start);

    // srand(1234);
    for (unsigned long i = 0; i < count; ++i) {
        // pairs_in.seek_relative((rand() % (pair_count / pair_storage_size)) *
        //                        pair_storage_size);

        WeightSystemPair pair{};

        read(pairs_in, pair.first);
        read(pairs_in, pair.second);

        weight_systems_from_pair(weight_systems, pair, statistics);

        // unsigned long new_candidate_count = candidate_count +
        // weight_systems.size();
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
        for (const auto &ws : weight_systems) {
            print_with_denominator(cout, ws);
            cout << endl;
        }
}

void check_ip(BufferedReader &ws_in)
{
    Stopwatch stopwatch{};

    check_config(ws_in);

    uint64_t ws_count;
    read(ws_in, ws_count);

    unsigned long ip_count = 0;
    for (unsigned long i = 0; i < ws_count; ++i) {
        WeightSystem<dim> ws{};
        read_varint(ws_in, ws);

        if (has_ip(ws)) {
            ++ip_count;
            if (g_settings.print_weight_systems) {
                print_with_denominator(cout, ws);
                cout << endl;
            }
        }
    }

    cerr << stopwatch << " - weight systems: " << ws_count
         << ", ip: " << ip_count << endl;
}

void read_ws(BufferedReader &ws_in)
{
    Stopwatch stopwatch{};

    check_config(ws_in);

    uint64_t ws_count;
    read(ws_in, ws_count);

    for (unsigned long i = 0; i < ws_count; ++i) {
        WeightSystem<dim> ws{};
        read_varint(ws_in, ws);

        if (g_settings.print_weight_systems) {
            print_with_denominator(cout, ws);
            cout << endl;
        }
    }

    cerr << stopwatch << " - weight systems: " << ws_count << endl;
}

void combine_ws_files(span<BufferedReader> ins, BufferedWriter &out)
{
    using std::pair;

    Stopwatch stopwatch{};

    vector<unsigned long> counts{};
    unsigned long total = 0;

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

    uint64_t count = total;
    write64u(out, count);
    count = 0;

    std::set<pair<WeightSystem<dim>, size_t>> weight_systems;

    for (size_t i = 0; i < counts.size(); ++i) {
        if (counts[i] > 0) {
            WeightSystem<dim> ws{};
            read_varint(ins[i], ws);
            weight_systems.insert(pair<WeightSystem<dim>, size_t>(ws, i));
        }
    }

    WeightSystem<dim> prev;

    while (!weight_systems.empty()) {
        auto smallest = *weight_systems.begin();
        weight_systems.erase(weight_systems.begin());

        const auto &smallest_ws = smallest.first;
        const auto &smallest_index = smallest.second;

        if (count == 0 || prev != smallest_ws) {
            ++count;
            write_varint(out, smallest_ws);
            prev = smallest_ws;
        }

        if (--counts[smallest_index] > 0) {
            WeightSystem<dim> ws{};
            read_varint(ins[smallest_index], ws);
            weight_systems.insert(
                pair<WeightSystem<dim>, size_t>(ws, smallest_index));
        }
    }

    if (total != count) {
        // Write header again, now with correct count
        out.seek(0);
        write_config(out);
        write64u(out, count);
    }

    cerr << stopwatch << " - combined weight systems: " << count << endl;
}

void split_ws_file(BufferedReader &in, span<BufferedWriter> outs)
{
    Stopwatch stopwatch{};

    check_config(in);

    uint64_t count;
    read(in, count);
    cerr << stopwatch << " - input weight systems: " << count << endl;

    unsigned long part_size = count / outs.size();
    unsigned long leftover = count % outs.size();
    unsigned long count_check = count;
    for (unsigned long i = 0; i < static_cast<unsigned long>(outs.size());
         ++i) {
        write_config(outs[i]);
        uint64_t out_count = i < leftover ? part_size + 1 : part_size;
        count_check -= out_count;
        write64u(outs[i], out_count);
    }
    assert(count_check == 0);

    while (count) {
        for (auto &out : outs) {
            WeightSystem<dim> ws;
            read_varint(in, ws);
            write_varint(out, ws);

            --count;
            if (count == 0)
                break;
        }
    }

    cerr << stopwatch << " - " << outs.size() << " files written\n";
}

WeightSystem<dim> add_half(const WeightSystem<dim - 1> &ws)
{
    assert(r_numerator == 1);
    assert(r_denominator == 1);

    WeightSystem<dim> ret{};
    for (unsigned j = 0; j < dim - 1; ++j)
        ret.weights[j] = ws.weights[j];

    ret.weights[dim - 1] = norm(ws);

    return ret;
}

void read_varint_lower(BufferedReader &in, WeightSystem<dim> &ws)
{
    WeightSystem<dim - 1> lower;
    read_varint(in, lower);
    ws = add_half(lower);
}

void split_two_ws_files(BufferedReader &in1, BufferedReader &in2,
                        span<BufferedWriter> outs)
{
    assert(r_denominator == 1);

    Stopwatch stopwatch{};

    check_config(in1);
    uint64_t count1;
    read(in1, count1);
    cerr << stopwatch << " - input weight systems: " << count1 << endl;

    check_config(in2, dim - 1, 2 * r_numerator - 1, 2);
    uint64_t count2;
    read(in2, count2);
    cerr << stopwatch << " - lower dimensional input weight systems: " << count2
         << endl;

    unsigned long total = count1 + count2;

    unsigned long part_size = total / outs.size();
    unsigned long leftover = total % outs.size();
    unsigned long count_check = total;
    for (unsigned long i = 0; i < static_cast<unsigned long>(outs.size());
         ++i) {
        write_config(outs[i]);
        uint64_t out_count = i < leftover ? part_size + 1 : part_size;
        count_check -= out_count;
        write64u(outs[i], out_count);
    }
    assert(count_check == 0);

    WeightSystem<dim> ws1{};
    WeightSystem<dim> ws2{};
    if (count1 > 0)
        read_varint(in1, ws1);
    if (count2 > 0)
        read_varint_lower(in2, ws2);

    while (count1 > 0 || count2 > 0) {
        for (auto &out : outs) {
            if ((count1 > 0 && count2 > 0 && ws1 < ws2) || count2 == 0) {
                write_varint(out, ws1);
                if (--count1 > 0)
                    read_varint(in1, ws1);
            } else {
                write_varint(out, ws2);
                if (--count2 > 0)
                    read_varint_lower(in2, ws2);
            }

            if (count1 == 0 && count2 == 0)
                break;
        }
    }

    cerr << stopwatch << " - " << outs.size() << " files written\n";
}

void diff_ws_files(BufferedReader &in1, BufferedReader &in2)
{
    Stopwatch stopwatch{};

    vector<unsigned long> counts{};

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

void analyze(BufferedReader &in, BufferedReader *info_in,
             BufferedWriter *info_out)
{
    Stopwatch stopwatch{};
    PolytopeStatistics stats{};

    check_config(in);

    uint64_t count;
    read(in, count);

    cerr << stopwatch << " - input weight systems: " << count << endl;

    for (unsigned long i = 0; i < count; ++i) {
        WeightSystem<dim> ws{};
        read_varint(in, ws);

        PolytopeInfo info{};

        if (info_in)
            read(*info_in, info);
        else
            analyze(ws, info, stats);

        if (info_out) {
            write(*info_out, info);
        } else {
            print_with_denominator(cout, ws);
            cout << " " << info << endl;
        }
    }

    if (info_in)
        cerr << stopwatch << endl;
    else
        cerr << stopwatch << " - statistics:\n" << stats;
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
    SwitchArg analyze_ws_arg( //
        "", "analyze-ws", "Analyze polytopes from given weight system file");
    SwitchArg split_ws_arg( //
        "", "split-ws", "Split weight system file into multiple files");
    SwitchArg read_ws_arg( //
        "", "read-ws", "Read weight systems from file");

    vector<Arg *> arg_list;
    arg_list.push_back(&find_candidates_arg);
    arg_list.push_back(&find_ip_arg);
    arg_list.push_back(&find_pairs_arg);
    arg_list.push_back(&print_count_arg);
    arg_list.push_back(&combine_ws_arg);
    arg_list.push_back(&diff_ws_arg);
    arg_list.push_back(&analyze_ws_arg);
    arg_list.push_back(&split_ws_arg);
    arg_list.push_back(&read_ws_arg);
    cmd.xorAdd(arg_list);

    SwitchArg print_ws_arg( //
        "", "print-ws", "Enable printing of weight systems", cmd);
    ValueArg<string> ws_out_arg( //
        "", "ws-out", "Weight systems destination file", false, "", "file",
        cmd);
    ValueArg<string> ws_in_arg( //
        "", "ws-in", "Weight systems source file", false, "", "file", cmd);
    ValueArg<string> ws_lower_in_arg( //
        "", "ws-lower-in", "Lower dimensional weight systems source file",
        false, "", "file", cmd);
    ValueArg<string> pairs_in_arg( //
        "", "pairs-in", "Pairs source file", false, "", "file", cmd);
    ValueArg<string> pairs_out_arg( //
        "", "pairs-out", "Pairs destination file", false, "", "file", cmd);
    ValueArg<unsigned long> start_arg( //
        "", "start", "", false, 0, "number", cmd);
    ValueArg<long> count_arg( //
        "", "count", "", false, -1, "number", cmd);
    ValueArg<string> polytope_info_out_arg( //
        "", "polytope-info-out", "Polytope info destination file", false, "",
        "file", cmd);
    ValueArg<string> polytope_info_in_arg( //
        "", "polytope-info-in", "Polytope info source file", false, "", "file",
        cmd);

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

    unsigned long start = start_arg.getValue();
    optional<unsigned long> count{};
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
    } else if (analyze_ws_arg.isSet() && ws_in_arg.isSet()) {
        BufferedReader in(ws_in_arg.getValue());

        unique_ptr<BufferedWriter> info_out{};
        if (polytope_info_out_arg.isSet())
            info_out =
                make_unique<BufferedWriter>(polytope_info_out_arg.getValue());

        unique_ptr<BufferedReader> info_in{};
        if (polytope_info_in_arg.isSet())
            info_in =
                make_unique<BufferedReader>(polytope_info_in_arg.getValue());

        analyze(in, info_in.get(), info_out.get());
    } else if (split_ws_arg.isSet() && ws_in_arg.isSet()) {
        BufferedReader in(ws_in_arg.getValue());

        vector<BufferedWriter> out{};
        for (const auto &filename : files_arg.getValue())
            out.push_back(BufferedWriter{filename});

        if (ws_lower_in_arg.isSet()) {
            BufferedReader lower_in{ws_lower_in_arg.getValue()};
            split_two_ws_files(in, lower_in, span<BufferedWriter>(out));
        } else {
            split_ws_file(in, span<BufferedWriter>(out));
        }
    } else if (read_ws_arg.isSet() && ws_in_arg.isSet()) {
        BufferedReader in(ws_in_arg.getValue());
        read_ws(in);
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
