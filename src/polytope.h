#ifndef POLYTOPE_H
#define POLYTOPE_H

#include "buffered_writer.h"
#include "vector.h"
#include "weight_system.h"
#include <array>
#include <limits>
#include <ostream>

struct PolytopeInfo {
    bool ip{};
    bool reflexive{};
    unsigned vertex_count{};
    unsigned facet_count{};
    unsigned dual_point_count{};
    unsigned point_count{};
    Vector<unsigned, dim - 3> hodge_numbers; // h1i
};

struct PolytopeStatistics {
    static constexpr int max_log2 = sizeof(Ring) * 8;

    unsigned long n_nonIP{};     // number of non IP weight systems
    unsigned long n_IP_nonRef{}; // number of non-reflexive IP weight systems
    unsigned long n_ref{};       // number of reflexive IP weight systems
    Ring max_w{};                // maximum weight in the reflexive cases
    Ring nr_max_w{};             // maximum weight in the non-reflexive cases
    std::array<unsigned long, max_log2>
        nr_n_w{}; // numbers of weights of given log2
    std::array<unsigned long, max_log2>
        n_w{};                         // numbers of weights of given log2
    std::array<int, dim - 1> max_nf{}; // max facet numbers

    // max values of certain entries of BH
    int nr_max_mp{};
    int nr_max_mv{};
    int nr_max_nv{};
    int max_mp{};
    int max_mv{};
    int max_np{};
    int max_nv{};
    int max_h22{};
    std::array<int, dim - 2> max_h1{};

    // range for chi (Euler number)
    int min_chi{std::numeric_limits<int>::max()};
    int max_chi{std::numeric_limits<int>::min()};
};

int euler_number(const PolytopeInfo &info);
int hodge_number_22(const PolytopeInfo &info);
void analyze(const WeightSystem<dim> &ws, PolytopeInfo &info,
             PolytopeStatistics &stats);
std::ostream &operator<<(std::ostream &os, const PolytopeInfo &info);
void read(BufferedReader &f, PolytopeInfo &info);
void write(BufferedWriter &f, const PolytopeInfo &info);

std::ostream &operator<<(std::ostream &os, const PolytopeStatistics &stats);

#endif
