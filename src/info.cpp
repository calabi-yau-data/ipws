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
using std::endl;
using std::make_unique;
using std::string;
using std::unique_ptr;
using std::unordered_set;
using std::vector;

struct WeightSystemWithInfo {
    WeightSystem<dim> ws{};
    PolytopeInfo info{};
};

bool operator<(const WeightSystemWithInfo &lhs, const WeightSystemWithInfo &rhs)
{
    if (lhs.info.ip != rhs.info.ip)
        return lhs.info.ip < rhs.info.ip;

    if (lhs.info.reflexive != rhs.info.reflexive)
        return lhs.info.reflexive < rhs.info.reflexive;

    if (lhs.info.reflexive) {
        if (lhs.info.hodge_numbers != rhs.info.hodge_numbers)
            return lhs.info.hodge_numbers < rhs.info.hodge_numbers;
    }

    if (lhs.info.vertex_count != rhs.info.vertex_count)
        return lhs.info.vertex_count < rhs.info.vertex_count;
    if (lhs.info.facet_count != rhs.info.facet_count)
        return lhs.info.facet_count < rhs.info.facet_count;
    if (lhs.info.point_count != rhs.info.point_count)
        return lhs.info.point_count < rhs.info.point_count;

    if (lhs.info.reflexive) {
        if (lhs.info.dual_point_count != rhs.info.dual_point_count)
            return lhs.info.dual_point_count < rhs.info.dual_point_count;
    }

    return lhs.ws < rhs.ws;
}

struct WeightSystemGroup {
    PolytopeInfo common_info{};
    std::vector<WeightSystemWithInfo> ws_list{};
};

void set_common_info(WeightSystemGroup &group, const PolytopeInfo &info)
{
    group.common_info.ip = info.ip;
    group.common_info.reflexive = info.reflexive;
    if (!info.ip)
        return;

    if (info.reflexive) {
        group.common_info.hodge_numbers = info.hodge_numbers;
    } else {
        group.common_info.vertex_count = info.vertex_count;
        group.common_info.facet_count = info.facet_count;
        group.common_info.point_count = info.point_count;
    }
}

bool same_common_info(const WeightSystemGroup &group, const PolytopeInfo &info)
{
    if (group.common_info.ip != info.ip)
        return false;

    if (!info.ip)
        return true;

    if (group.common_info.reflexive != info.reflexive)
        return false;

    if (info.reflexive) {
        return group.common_info.hodge_numbers == info.hodge_numbers;
    } else {
        return group.common_info.vertex_count == info.vertex_count &&
               group.common_info.facet_count == info.facet_count &&
               group.common_info.point_count == info.point_count;
    }
}

bool same_common_info(const WeightSystemGroup &group,
                      const WeightSystemGroup &group2)
{
    return same_common_info(group, group2.common_info);
}

bool operator<(const WeightSystemGroup &lhs, const WeightSystemGroup &rhs)
{
    if (lhs.common_info.ip != rhs.common_info.ip)
        return lhs.common_info.ip < rhs.common_info.ip;

    if (!lhs.common_info.ip)
        return false;

    if (lhs.common_info.reflexive != rhs.common_info.reflexive)
        return lhs.common_info.reflexive < rhs.common_info.reflexive;

    if (lhs.common_info.reflexive) {
        return lhs.common_info.hodge_numbers < rhs.common_info.hodge_numbers;
    } else {
        if (lhs.common_info.vertex_count != rhs.common_info.vertex_count)
            return lhs.common_info.vertex_count < rhs.common_info.vertex_count;
        if (lhs.common_info.facet_count != rhs.common_info.facet_count)
            return lhs.common_info.facet_count < rhs.common_info.facet_count;
        return lhs.common_info.point_count < rhs.common_info.point_count;
    }
}

void read(BufferedReader &f, WeightSystemGroup &group, bool reflexive)
{
    group.common_info.ip = true;
    group.common_info.reflexive = reflexive;

    if (reflexive) {
        read_varint(f, group.common_info.hodge_numbers);
        group.ws_list.resize(read_varint(f));

        for (auto &x : group.ws_list) {
            read_varint(f, x.ws);
            x.info.vertex_count = static_cast<unsigned>(read_varint(f));
            x.info.facet_count = static_cast<unsigned>(read_varint(f));
            x.info.point_count = static_cast<unsigned>(read_varint(f));
            x.info.dual_point_count = static_cast<unsigned>(read_varint(f));
        }
    } else {
        group.common_info.vertex_count = static_cast<unsigned>(read_varint(f));
        group.common_info.facet_count = static_cast<unsigned>(read_varint(f));
        group.common_info.point_count = static_cast<unsigned>(read_varint(f));
        group.ws_list.resize(read_varint(f));

        for (auto &x : group.ws_list)
            read_varint(f, x.ws);
    }
}

void write(BufferedWriter &f, const WeightSystemGroup &group)
{
    if (!group.common_info.ip)
        return;

    if (group.common_info.reflexive) {
        write_varint(f, group.common_info.hodge_numbers);
        write_varint(f, group.ws_list.size());

        for (const auto &x : group.ws_list) {
            write_varint(f, x.ws);
            write_varint(f, x.info.vertex_count);
            write_varint(f, x.info.facet_count);
            write_varint(f, x.info.point_count);
            write_varint(f, x.info.dual_point_count);
        }
    } else {
        write_varint(f, group.common_info.vertex_count);
        write_varint(f, group.common_info.facet_count);
        write_varint(f, group.common_info.point_count);
        write_varint(f, group.ws_list.size());

        for (const auto &x : group.ws_list)
            write_varint(f, x.ws);
    }
}

void write_grouped(BufferedWriter &out,
                   const vector<WeightSystemWithInfo> &list)
{
    write_config(out);

    for (auto it = list.begin(); it != list.end();) {
        WeightSystemGroup group{};
        set_common_info(group, it->info);

        auto block_end = it;
        while (block_end != list.end() &&
               same_common_info(group, block_end->info))
            ++block_end;

        group.ws_list.resize(block_end - it);
        std::copy(it, block_end, group.ws_list.begin());

        write(out, group);

        it = block_end;
    }
}

void sort_with_info(BufferedReader &ws_in, BufferedReader &info_in,
                    BufferedWriter *reflexive_out,
                    BufferedWriter *non_reflexive_out)
{
    Stopwatch stopwatch{};

    check_config(ws_in);

    uint64_t count;
    read(ws_in, count);

    cerr << stopwatch << " - input weight systems: " << count << endl;
    if (count == 0)
        return;

    vector<WeightSystemWithInfo> reflexive_list{};
    vector<WeightSystemWithInfo> non_reflexive_list{};

    for (unsigned long i = 0; i < count; ++i) {
        WeightSystemWithInfo x{};
        read_varint(ws_in, x.ws);
        read(info_in, x.info);

        if (!x.info.ip)
            continue;

        if (x.info.reflexive && reflexive_out)
            reflexive_list.push_back(x);

        if (!x.info.reflexive && non_reflexive_out)
            non_reflexive_list.push_back(x);
    }

    cerr << stopwatch << " - done reading\n";

    std::sort(reflexive_list.begin(), reflexive_list.end());
    std::sort(non_reflexive_list.begin(), non_reflexive_list.end());

    cerr << stopwatch << " - done sorting\n";

    if (reflexive_out) {
        write_grouped(*reflexive_out, reflexive_list);
        cerr << stopwatch << " - reflexive count: " << reflexive_list.size()
             << endl;
    }

    if (non_reflexive_out) {
        write_grouped(*non_reflexive_out, non_reflexive_list);
        cerr << stopwatch
             << " - non-reflexive count: " << non_reflexive_list.size() << endl;
    }
}

void combine(BufferedWriter &out, span<BufferedReader> ins, bool reflexive)
{
    using std::pair;

    Stopwatch stopwatch{};

    std::set<pair<WeightSystemGroup, size_t>> group_list{};

    for (size_t i = 0; i < static_cast<size_t>(ins.size()); ++i) {
        check_config(ins[i]);
        try {
            WeightSystemGroup g{};
            read(ins[i], g, reflexive);
            group_list.insert(pair<WeightSystemGroup, size_t>(g, i));
        } catch (BufferedReader::EofError) {
        }
    }

    WeightSystemGroup prev{};
    bool prev_set = false;
    unsigned long count = 0;
    unsigned long group_count = 0;

    write_config(out);

    while (!group_list.empty()) {
        auto smallest = *group_list.begin();
        group_list.erase(group_list.begin());

        const auto &smallest_group = smallest.first;
        const auto &smallest_index = smallest.second;

        if (prev_set) {
            bool max_size_reached =
                prev.ws_list.size() + smallest_group.ws_list.size() > 1000000;

            if (!max_size_reached && same_common_info(prev, smallest_group)) {
                std::copy(smallest_group.ws_list.begin(),
                          smallest_group.ws_list.end(),
                          std::back_inserter(prev.ws_list));
            } else {
                ++group_count;
                count += prev.ws_list.size();

                write(out, prev);

                prev = smallest_group;

                if (group_count % 1000000 == 0) {
                    cerr << stopwatch << " - weight system count: " << count
                         << endl;
                    cerr << stopwatch << " - groups: " << group_count << endl;
                }
            }
        } else {
            prev = smallest_group;
            prev_set = true;
        }

        try {
            WeightSystemGroup g{};
            read(ins[smallest_index], g, reflexive);
            group_list.insert(
                pair<WeightSystemGroup, size_t>(g, smallest_index));
        } catch (BufferedReader::EofError) {
        }
    }

    if (prev_set) {
        ++group_count;
        count += prev.ws_list.size();

        write(out, prev);
    }

    cerr << stopwatch << " - weight system count: " << count << endl;
    cerr << stopwatch << " - groups: " << group_count << endl;
}

/*
https://www.postgresql.org/docs/9.6/static/sql-copy.html
create table reflexive (
  h11 int4 not null,
  h12 int4 not null,
  h13 int4 not null,
  h22 int4 not null,
  chi int4 not null,
  ws_count int4 not null,
  ws_data bytea not null
);
copy reflexive from 'pgcopy-file-path' binary;
*/
void pgcopy(BufferedReader &in, BufferedWriter &out, bool reflexive)
{
    assert(reflexive);
    assert(dim == 6);

    Stopwatch stopwatch{};

    check_config(in);

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

    unsigned long group_count = 0;
    unsigned long ws_count = 0;

    for (;;) {
        WeightSystemGroup g{};
        try {
            read(in, g, reflexive);
        } catch (BufferedReader::EofError) {
            break;
        }

        ++group_count;
        ws_count += g.ws_list.size();

        write16i(out,
                 static_cast<int16_t>(g.common_info.hodge_numbers.size() + 4));

        for (const auto &n : g.common_info.hodge_numbers.array) {
            write32i(out, 4);
            write32i(out, n);
        }
        write32i(out, 4);
        write32i(out, picard_number(g.common_info));
        write32i(out, 4);
        write32i(out, euler_number(g.common_info));
        write32i(out, 4);
        write32i(out, static_cast<int32_t>(g.ws_list.size()));

        unsigned size = 0;
        for (const auto &x : g.ws_list) {
            for (const auto &weight : x.ws.weights)
                size += varint_storage_size(weight);
            size += varint_storage_size(x.info.vertex_count);
            size += varint_storage_size(x.info.facet_count);
            size += varint_storage_size(x.info.point_count);
            size += varint_storage_size(x.info.dual_point_count);
        }

        write32i(out, size);
        for (const auto &x : g.ws_list) {
            for (const auto &weight : x.ws.weights)
                write_varint(out, weight);
            write_varint(out, x.info.vertex_count);
            write_varint(out, x.info.facet_count);
            write_varint(out, x.info.point_count);
            write_varint(out, x.info.dual_point_count);
        }
    }

    write16i(out, -1);

    cerr << stopwatch << " - group count: " << group_count << endl;
    cerr << stopwatch << " - weight systems: " << ws_count << endl;
}

bool run(int argc, char *argv[])
{
    using TCLAP::Arg;
    using TCLAP::SwitchArg;
    using TCLAP::ValueArg;
    using std::fstream;

    std::ostringstream description;
    description << "Classify weight systems with d=" << dim
                << ", r=" << r_numerator << "/" << r_denominator
                << " using ring '" << ring_name << "'";

    TCLAP::CmdLine cmd(description.str(), ' ', GIT_REVISION);

    SwitchArg sort_arg( //
        "", "sort", "Sort according to polytope info");
    SwitchArg combine_arg( //
        "", "combine", "");
    ValueArg<string> pgcopy_arg( //
        "", "pgcopy", "PostgreSQL COPY destination file", false, "", "file");

    vector<Arg *> arg_list;
    arg_list.push_back(&sort_arg);
    arg_list.push_back(&combine_arg);
    arg_list.push_back(&pgcopy_arg);
    cmd.xorAdd(arg_list);

    ValueArg<string> ws_in_arg( //
        "", "ws-in", "Weight systems source file", false, "", "file", cmd);
    ValueArg<string> polytope_info_in_arg( //
        "", "polytope-info-in", "Polytope info source file", false, "", "file",
        cmd);

    ValueArg<string> reflexive_in_arg( //
        "", "reflexive-in", "Reflexive source file", false, "", "file", cmd);
    ValueArg<string> non_reflexive_in_arg( //
        "", "non-reflexive-in", "Non-reflexive source file", false, "", "file",
        cmd);
    ValueArg<string> reflexive_out_arg( //
        "", "reflexive-out", "Reflexive destination file", false, "", "file",
        cmd);
    ValueArg<string> non_reflexive_out_arg( //
        "", "non-reflexive-out", "Non-reflexive destination file", false, "",
        "file", cmd);

    TCLAP::UnlabeledMultiArg<string> files_arg( //
        "", "Additional files", false, "file", cmd);

    cmd.parse(argc, argv);

    if (sort_arg.isSet() && ws_in_arg.isSet() && polytope_info_in_arg.isSet()) {
        unique_ptr<BufferedWriter> reflexive_out{};
        if (reflexive_out_arg.isSet())
            reflexive_out =
                make_unique<BufferedWriter>(reflexive_out_arg.getValue());

        unique_ptr<BufferedWriter> non_reflexive_out{};
        if (non_reflexive_out_arg.isSet())
            non_reflexive_out =
                make_unique<BufferedWriter>(non_reflexive_out_arg.getValue());

        BufferedReader ws_in(ws_in_arg.getValue());
        BufferedReader info_in(polytope_info_in_arg.getValue());

        sort_with_info(ws_in, info_in, reflexive_out.get(),
                       non_reflexive_out.get());
    } else if (combine_arg.isSet()) {
        vector<BufferedReader> ins{};
        for (const auto &filename : files_arg.getValue())
            ins.push_back(BufferedReader(filename));

        if (reflexive_out_arg.isSet() && !non_reflexive_out_arg.isSet()) {
            BufferedWriter out(reflexive_out_arg.getValue());
            combine(out, ins, true);
        }

        if (!reflexive_out_arg.isSet() && non_reflexive_out_arg.isSet()) {
            BufferedWriter out(non_reflexive_out_arg.getValue());
            combine(out, ins, false);
        }
    } else if (pgcopy_arg.isSet()) {
        BufferedWriter out(pgcopy_arg.getValue());

        if (reflexive_in_arg.isSet() && !non_reflexive_in_arg.isSet()) {
            BufferedReader in(reflexive_in_arg.getValue());
            pgcopy(in, out, true);
        }

        if (!reflexive_in_arg.isSet() && non_reflexive_in_arg.isSet()) {
            BufferedReader in(non_reflexive_in_arg.getValue());
            pgcopy(in, out, false);
        }
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
