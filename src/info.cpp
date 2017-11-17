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
#include "polytope.h"
#include "read_write.h"
#include "settings.h"
#include "stl_utils.h"
#include "stopwatch.h"
#include "weight_system.h"
#include "weight_system_builder.h"
#include "weight_system_pair.h"

using boost::optional;
using gsl::span;
using std::array;
using std::cerr;
using std::cout;
using std::endl;
using std::unordered_set;
using std::string;
using std::vector;
using std::unique_ptr;
using std::make_unique;

struct HodgeNumbersWithCount {
    Vector<unsigned, dim - 3> hodge_numbers{};
    unsigned long count;

    bool operator<(const HodgeNumbersWithCount &rhs) const
    {
        for (unsigned i = 0; i < hodge_numbers.size(); ++i)
            if (hodge_numbers[i] != rhs.hodge_numbers[i])
                return hodge_numbers[i] < rhs.hodge_numbers[i];
        return count < rhs.count;
    }
};

void read_varint(BufferedReader &f, HodgeNumbersWithCount &x)
{
    read_varint(f, x.hodge_numbers);
    x.count = read_varint(f);
}

void write_varint(BufferedWriter &f, const HodgeNumbersWithCount &x)
{
    write_varint(f, x.hodge_numbers);
    write_varint(f, x.count);
}

void sql(BufferedReader &in, BufferedWriter &out)
{
    // https://www.postgresql.org/docs/9.6/static/sql-copy.html

    Stopwatch stopwatch{};

    write8u(out, 'P');
    write8u(out, 'G');
    write8u(out, 'C');
    write8u(out, 'O');
    write8u(out, 'P');
    write8u(out, 'Y');
    write8u(out, '\n');
    write8u(out, '\377');
    write8u(out, '\r');
    write8u(out, '\n');
    write8u(out, '\0');

    write32u(out, 0);
    write32i(out, 0);

    unsigned long hodge_count = 0;
    unsigned long ws_count = 0;

    for (;;) {
        HodgeNumbersWithCount x{};

        try {
            read_varint(in, x);
        } catch (BufferedReader::EofError) {
            break;
        }

        ++hodge_count;
        ws_count += x.count;

        write16i(out, static_cast<int16_t>(x.hodge_numbers.size() + 1));

        for (const auto &n : x.hodge_numbers.array) {
            write32i(out, 4);
            write32i(out, n);
        }
        write32i(out, 4);
        write32i(out, static_cast<int32_t>(x.count));
    }

    write16i(out, -1);

    cerr << stopwatch << " - hodge numbers: " << hodge_count << endl;
    cerr << stopwatch << " - weight systems: " << ws_count << endl;
}

struct WeightSystemWithHodgeNumbers {
    WeightSystem<dim> ws{};
    Vector<unsigned, dim - 3> hodge_numbers{};
};

bool operator<(const WeightSystemWithHodgeNumbers &lhs,
               const WeightSystemWithHodgeNumbers &rhs)
{
    for (unsigned i = 0; i < lhs.hodge_numbers.size(); ++i)
        if (lhs.hodge_numbers[i] != rhs.hodge_numbers[i])
            return lhs.hodge_numbers[i] < rhs.hodge_numbers[i];

    return compare(lhs.ws, rhs.ws) < 0;
}

void sort_hodge(BufferedReader &ws_in, BufferedReader &info_in,
                BufferedWriter &ws_out, BufferedWriter &info_out)
{
    Stopwatch stopwatch{};

    check_config(ws_in);

    uint64_t count;
    read(ws_in, count);

    cerr << stopwatch << " - input weight systems: " << count << endl;
    if (count == 0)
        return;

    vector<WeightSystemWithHodgeNumbers> ws_list{};
    ws_list.reserve(count);

    for (unsigned long i = 0; i < count; ++i) {
        WeightSystemWithHodgeNumbers x{};
        read_varint(ws_in, x.ws);

        PolytopeInfo info{};
        read(info_in, info);

        if (!info.ip || !info.reflexive)
            continue;

        x.hodge_numbers = info.hodge_numbers;

        ws_list.push_back(x);
    }

    cerr << stopwatch << " - done reading\n";

    std::sort(ws_list.begin(), ws_list.end());

    cerr << stopwatch << " - done sorting\n";

    write_config(ws_out);
    write64u(ws_out, ws_list.size());

    unsigned long same_hodge_count = 0;
    for (unsigned int i = 0; i < ws_list.size(); ++i) {
        const auto &current = ws_list[i];
        if (i > 0) {
            const auto &prev = ws_list[i - 1];
            if (prev.hodge_numbers != current.hodge_numbers) {
                write_varint(info_out, prev.hodge_numbers);
                write_varint(info_out, same_hodge_count);
                same_hodge_count = 0;
            }
        }

        ++same_hodge_count;
        write_varint(ws_out, current.ws);
    }

    write_varint(info_out, ws_list.back().hodge_numbers);
    write_varint(info_out, same_hodge_count);

    cerr << stopwatch << endl;
}

// TODO: check this function!
void combine_hodge(BufferedWriter &info_out, span<BufferedReader> info_ins)
{
    using std::pair;

    Stopwatch stopwatch{};

    std::set<pair<HodgeNumbersWithCount, size_t>> hodge_numbers_list{};

    for (size_t i = 0; i < static_cast<size_t>(info_ins.size()); ++i) {
        try {
            HodgeNumbersWithCount x{};
            read_varint(info_ins[i], x);
            hodge_numbers_list.insert(
                pair<HodgeNumbersWithCount, size_t>(x, i));
        } catch (BufferedReader::EofError) {
        }
    }

    HodgeNumbersWithCount prev;
    unsigned long count = 0;
    unsigned long unique_count = 0;
    bool first_iteration = true;

    while (!hodge_numbers_list.empty()) {
        auto smallest = *hodge_numbers_list.begin();
        hodge_numbers_list.erase(hodge_numbers_list.begin());

        const auto &smallest_hodge = smallest.first;
        const auto &smallest_index = smallest.second;

        if (!first_iteration) {
            if (prev.hodge_numbers == smallest_hodge.hodge_numbers) {
                prev.count += smallest_hodge.count;
            } else {
                ++unique_count;
                count += prev.count;

                write_varint(info_out, prev);

                prev = smallest_hodge;

                if (unique_count % 1000000 == 0) {
                    cerr << stopwatch << " - weight system count: " << count
                         << endl;
                    cerr << stopwatch << " - unique hodge number combinations: "
                         << unique_count << endl;
                }
            }
        }
        first_iteration = false;

        try {
            HodgeNumbersWithCount x{};
            read_varint(info_ins[smallest_index], x);
            hodge_numbers_list.insert(
                pair<HodgeNumbersWithCount, size_t>(x, smallest_index));
        } catch (BufferedReader::EofError) {
        }
    }

    if (!first_iteration) {
        ++unique_count;
        count += prev.count;

        write_varint(info_out, prev);
    }

    cerr << stopwatch << " - weight system count: " << count << endl;
    cerr << stopwatch << " - unique hodge number combinations: " << unique_count
         << endl;
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

    SwitchArg sort_hodge_arg( //
        "", "sort-hodge", "Sort according to hodge numbers");
    SwitchArg combine_hodge_arg( //
        "", "combine-hodge", "");
    SwitchArg sql_arg( //
        "", "sql", "");

    vector<Arg *> arg_list;
    arg_list.push_back(&sort_hodge_arg);
    arg_list.push_back(&combine_hodge_arg);
    arg_list.push_back(&sql_arg);
    cmd.xorAdd(arg_list);

    ValueArg<string> ws_out_arg( //
        "", "ws-out", "Weight systems destination file", false, "", "file",
        cmd);
    ValueArg<string> ws_in_arg( //
        "", "ws-in", "Weight systems source file", false, "", "file", cmd);
    ValueArg<string> polytope_info_out_arg( //
        "", "polytope-info-out", "Polytope info destination file", false, "",
        "file", cmd);
    ValueArg<string> polytope_info_in_arg( //
        "", "polytope-info-in", "Polytope info source file", false, "", "file",
        cmd);

    TCLAP::UnlabeledMultiArg<string> files_arg( //
        "", "Additional files", false, "file", cmd);

    cmd.parse(argc, argv);

    if (sort_hodge_arg.isSet() && ws_in_arg.isSet() &&
        polytope_info_in_arg.isSet() && ws_out_arg.isSet() &&
        polytope_info_out_arg.isSet()) {
        BufferedReader ws_in(ws_in_arg.getValue());
        BufferedWriter ws_out(ws_out_arg.getValue());
        BufferedReader info_in(polytope_info_in_arg.getValue());
        BufferedWriter info_out(polytope_info_out_arg.getValue());
        sort_hodge(ws_in, info_in, ws_out, info_out);
    } else if (combine_hodge_arg.isSet() && polytope_info_out_arg.isSet()) {
        vector<BufferedReader> info_ins{};
        for (const auto &filename : files_arg.getValue())
            info_ins.push_back(BufferedReader(filename));

        BufferedWriter info_out(polytope_info_out_arg.getValue());

        combine_hodge(info_out, info_ins);
    } else if (sql_arg.isSet() && polytope_info_in_arg.isSet() &&
               polytope_info_out_arg.isSet()) {
        BufferedReader in(polytope_info_in_arg.getValue());
        BufferedWriter out(polytope_info_out_arg.getValue());

        sql(in, out);
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
