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

void write_config(BufferedWriter &f, unsigned dim_, Ring r_numerator_,
                  Ring r_denominator_)
{
    write(f, static_cast<uint32_t>(dim_));
    write(f, static_cast<uint32_t>(r_numerator_));
    write(f, static_cast<uint32_t>(r_denominator_));
}

void add_half(BufferedReader &in, BufferedWriter &out)
{
    assert(r_numerator % 2 == 1);
    assert(r_denominator == 2);

    check_config(in);
    uint64_t count;
    read(in, count);

    write_config(out, dim + 1, (r_numerator + 1) / 2, 1);
    write(out, count);

    for (size_t i = 0; i < count; ++i) {
        WeightSystem<dim> ws_in;
        WeightSystem<dim + 1> ws_out;

        read_varint(in, ws_in);

        for (unsigned j = 0; j < dim; ++j)
            ws_out.weights[j] = ws_in.weights[j];

        assert(norm(ws_in) % r_numerator == 0);
        ws_out.weights[dim] = norm(ws_in) / r_numerator;

        write_varint(out, ws_out);

        // cout << ws_in << " -> " << ws_out << endl;
    }
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

    if (ws_in_arg.isSet() && ws_out_arg.isSet()) {
        BufferedReader ws_in(ws_in_arg.getValue());
        BufferedWriter ws_out(ws_out_arg.getValue());
        add_half(ws_in, ws_out);
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
