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

constexpr unsigned every = 10000;

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

void write_config(BufferedWriter &f)
{
    write(f, static_cast<uint32_t>(dim));
    write(f, static_cast<uint32_t>(r_numerator));
    write(f, static_cast<uint32_t>(r_denominator));
}

void print_with_denominator(const WeightSystem<dim> &ws)
{
    Ring n = norm(ws);
    cout << n * r_denominator;
    for (const auto &w : ws.weights)
        cout << " " << w * r_numerator;
    cout << endl;
}

void sample(BufferedReader &in)
{
    check_config(in);
    uint64_t count;
    read(in, count);

    for (unsigned long i = 0; i < count; ++i) {
        WeightSystem<dim> ws;

        read_varint(in, ws);
        if (i % every == 0)
            print_with_denominator(ws);
    }
}

void sample(BufferedReader &in, BufferedWriter &out)
{
    check_config(in);
    uint64_t count;
    read(in, count);

    write_config(out);
    uint64_t count_out = count / every + 1;
    write(out, count_out);

    unsigned long written = 0;
    for (unsigned long i = 0; i < count; ++i) {
        WeightSystem<dim> ws;

        read_varint(in, ws);
        if (i % every == 0) {
            write_varint(out, ws);
            ++written;
        }
    }
    assert(written == count_out);
}

bool run(int argc, char *argv[])
{
    using std::fstream;
    using TCLAP::Arg;
    using TCLAP::SwitchArg;
    using TCLAP::ValueArg;

    TCLAP::CmdLine cmd("", ' ', GIT_REVISION);

    vector<Arg *> arg_list;

    ValueArg<string> ws_in_arg( //
        "", "ws-in", "Weight systems source file", false, "", "file", cmd);
    ValueArg<string> ws_out_arg( //
        "", "ws-out", "Weight systems destination file", false, "", "file",
        cmd);

    cmd.parse(argc, argv);

    if (ws_in_arg.isSet()) {
        BufferedReader ws_in(ws_in_arg.getValue());

        if (ws_out_arg.isSet()) {
            BufferedWriter ws_out(ws_out_arg.getValue());
            sample(ws_in, ws_out);
        } else {
            sample(ws_in);
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
