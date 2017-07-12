#include <iostream>
#include <set>
#include "config.h"
#include "point.h"
#include "stl_utils.h"
#include "stopwatch.h"
#include "weight_system.h"
#include "weight_system_builder.h"
#include "weight_system_pair.h"

int Find_Equations(PolyPointList *P, VertexNumList *VNL, EqList *EL);
Long Eval_Eq_on_V(Equation *E, Long *V, int i);

using gsl::span;
using std::array;
using std::cerr;
using std::cout;
using std::endl;
using std::pair;
using std::set;
using std::setprecision;
using std::vector;

set<WeightSystemPair> pairs; // TODO: unordered_set? TODO: not global

struct History {
    array<Point, dim - 1> points;
    array<WeightSystem, dim> weight_systems;
    array<array<Long, dim - 1>, dim - 1> point_weight_system_distances;
};

struct Statistics {
    unsigned weight_systems_found;
    unsigned final_pairs_found;
    unsigned ip_weight_systems;
};

using WeightSystemCollection = set<WeightSystem>; // TODO: unordered_set?

bool good_weight_system(const WeightSystem &ws)
{
    Long n = norm(ws);

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
            Long diff = history.point_weight_system_distances[i][i] -
                        distance(history.weight_systems[i], x);
            if (diff > 0 ||
                (!disable_lex_compare && diff == 0 && x > history.points[i]))
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
        Long rel = history.point_weight_system_distances[i][i] -
                   history.point_weight_system_distances[n][i];
        if (rel > 0)
            return true;
        if (!disable_lex_compare && rel == 0 && history.points[i] < x)
            return true;
    }

    for (int i = 0; i < n; ++i) {
        Point x_other = history.points[i];
        Point x_diff = x_other - x;

        Long v = gcd(x_diff);
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

void rec(WeightSystemCollection &weight_systems,
         const WeightSystemBuilder &builder, unsigned n, History &history,
         Statistics &statistics)
{
    WeightSystem ws{};
    if (!builder.average_if_nonzero(ws))
        return;

    history.weight_systems[n] = ws;

    if (n < dim - redundancy_check_skip_recursions &&
        last_point_redundant2(builder, n, history))
        return;

    add_maybe(weight_systems, ws, statistics);

    switch (n) {
    case dim - 2:
        // The following happens when redundancies are not checked
        assert(builder.generator_count() == 2);
        // TODO: special case for last iteration brings a big performance boost
        if (defer_last_recursion) {
            pairs.insert(canonicalize(builder.to_pair()));
            ++statistics.final_pairs_found;

            // if (final_pairs_found % 10000 == 0)
            //     System.out.printf("%7.2f: builders: %d, unique: %d\n",
            //                       stopwatch.count(), buildersAddedCount,
            //                       builders.size());
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

        rec(weight_systems, builder.restrict(x), n + 1, history, statistics);
    }
}

void add_point(span<const Long> x, PolyPointList *P)
{
    assert(P->np < POINT_Nmax);
    for (unsigned i = 0; i < dim; i++)
        P->x[P->np][i] = x[i];
    P->np++;
}

PolyPointList *new_point_list(const WeightSystem &ws)
{
    PolyPointList *P = (PolyPointList *)malloc(sizeof(PolyPointList));
    assert(P != nullptr);

    P->np = 0;

    add_point(Point().coords, P);

    auto gen = WeightSystemPointsOn(ws);
    while (gen.find_next()) {
        add_point(gen.get().coords, P);
    }

    return P;
}

bool has_ip(const WeightSystem &ws)
{
    // static auto poly_point_list = std::make_unique<PolyPointList>();
    PolyPointList *P = new_point_list(ws);

    VertexNumList V;
    EqList E;

    if (static_cast<unsigned>(P->np) <= dim) {
        free(P);
        return 0;
    }
    Find_Equations(P, &V, &E);
    if (static_cast<unsigned>(E.ne) < dim) {
        free(P);
        return 0;
    }
    Long y[dim];
    for (unsigned k = 0; k < dim; k++)
        y[k] = 1;
    for (unsigned k = 0; k < static_cast<unsigned>(E.ne); k++)
        if (Eval_Eq_on_V(&(E.e[k]), y, dim) <= 0)
            if (!E.e[k].c) {
                free(P);
                return 0;
            }

    bool ret = P->np != 1;
    free(P);
    return ret;
}

int main()
{
    Stopwatch stopwatch{};

    History history{};
    WeightSystemCollection weight_systems{};
    Statistics statistics{};

    rec(weight_systems, WeightSystemBuilder{}, 0, history, statistics);

    if (defer_last_recursion) {
        cerr << stopwatch
             << " - weight systems: " << statistics.weight_systems_found
             << ", unique: " << weight_systems.size() << endl;
        cerr << stopwatch << " - pairs: " << statistics.final_pairs_found
             << ", unique: " << pairs.size() << endl;

        for (const auto &pair : pairs) {
            WeightSystem ws = average(pair);
            unsigned candidate_count = 0;
            unsigned old_unique_count = weight_systems.size();
            set<WeightSystem> current_weight_systems{}; // TODO: unordered_set?

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
                    if (add_maybe(weight_systems, final_ws, statistics) &&
                        print_last_recursion_statistics) {
                        ++candidate_count;
                        current_weight_systems.insert(final_ws);
                    }
                }
            }

            if (print_last_recursion_statistics) {
                unsigned unique_count =
                    weight_systems.size() - old_unique_count;
                cout << candidate_count << " " << current_weight_systems.size()
                     << " " << unique_count << endl;
            }
        }
    }

    cerr << stopwatch
         << " - weight systems: " << statistics.weight_systems_found
         << ", unique: " << weight_systems.size() << endl;

    for (const auto &ws : weight_systems) {
        if (has_ip(ws))
            ++statistics.ip_weight_systems;
    }

    cerr << stopwatch
         << " - weight systems: " << statistics.weight_systems_found
         << ", unique: " << weight_systems.size()
         << ", ip: " << statistics.ip_weight_systems << endl;

    return 0;
}
