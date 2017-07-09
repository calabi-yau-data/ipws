#include <arpa/inet.h>
#include <time.h>
#include <algorithm>
#include <array>
#include <bitset>
#include <chrono>
#include <climits>
#include <fstream>
#include <iostream>
#include <memory>
#include <set>
#include <vector>

#include "Global.h"
#include "Rat.h"
#include "hyperplane.h"
#include "vector.h"

#include "weight_system_store.h"

using std::cout;
using std::endl;

namespace {


const bool write_cones = false;
const bool read_from_file = false;
const bool defer_last_recursion = false;
const bool disable_lex_cmp = true;

// it is not economical to do this on the last two recursion levels
const int redundancy_check_skip_recursions = 2;

template <class T, class F>
void rearranging_erase_if(std::vector<T> &c, F f)
{
    for (size_t i = c.size(); i--;) {
        if (f(c[i])) {
            c[i] = c.back();
            c.pop_back();
        }
    }
}

template <size_t N>
bool operator<(const std::bitset<N> &lhs, const std::bitset<N> &rhs)
{
    return (lhs | rhs) == rhs;
}

// Enumerates all points x with nonnegative coordinates, satisfying
// q.a * x + q.c < 0. Very optimized.
template <class F>
void enumerate_points_below(const Hyperplane &q, F callback)
{
    std::array<Long, dim> x;
    std::array<Long, dim> ax;

    std::fill_n(x.begin(), dim - 1, 0);
    std::fill_n(ax.begin(), dim - 1, 0);

    x[dim - 1] = -1;
    ax[dim - 1] = -q.a[dim - 1];

    for (int k = dim - 1; k >= 0;) {
        x[k]++;
        ax[k] += q.a[k];
        for (int i = k + 1; i < dim; i++)
            ax[i] = ax[k];

        callback(x);

        for (k = dim - 1; k >= 0 && ax[k] + q.a[k] + q.c >= 0; k--)
            x[k] = 0;
    }
}

// Enumerates all points x with nonnegative coordinates, satisfying
// q.a * x + q.c == 0. Not optimized.
template <class F>
void enumerate_points_on(const Hyperplane &q, F callback)
{
    std::array<Long, dim> x;
    std::array<Long, dim> ax;

    std::fill_n(x.begin(), dim - 1, 0);
    std::fill_n(ax.begin(), dim - 1, 0);

    x[dim - 1] = -1;
    ax[dim - 1] = -q.a[dim - 1];

    for (int k = dim - 1; k >= 0;) {
        x[k]++;
        ax[k] += q.a[k];
        for (int i = k + 1; i < dim; i++)
            ax[i] = ax[k];

        if (ax[k] + q.c == 0)
            callback(x);

        for (k = dim - 1; k >= 0 && ax[k] + q.a[k] + q.c > 0; k--)
            x[k] = 0;
    }
}

struct RationalCone {
    struct Generator {
        Hyperplane eq;
        std::bitset<dim> incidences;
    };
    std::vector<Generator> generators;

public:
    bool operator<(const RationalCone &rhs) const
    {
        if (generators.size() != rhs.generators.size())
            return generators.size() < rhs.generators.size();

        for (size_t i = 0; i < generators.size(); ++i)
            if (generators[i].eq != rhs.generators[i].eq)
                return generators[i].eq < rhs.generators[i].eq;

        return false;
    }

    // c has to be smaller than 0 for this to really yield the positive orthant.
    static RationalCone positive_orthant(Long a, Long c)
    {
        RationalCone ret{};

        ret.generators.resize(dim);
        for (size_t generator_nr = 0; generator_nr < dim; ++generator_nr) {
            auto &generator = ret.generators[generator_nr];

            for (size_t i = 0; i < dim; ++i)
                generator.eq.a[i] = generator_nr == i ? a : 0;
            generator.eq.c = c;

            generator.incidences.set(generator_nr);
        }

        return ret;
    }

    RationalCone restrict(const Vector &x, int n) const
    {
        RationalCone ret{};

        std::vector<Long> xq;
        xq.resize(generators.size());

        for (unsigned int i = 0; i < generators.size(); i++) {
            xq[i] = generators[i].eq * x;
            if (xq[i] == 0)
                ret.generators.push_back(generators[i]);
        }

        for (unsigned int i = 0; i < generators.size() - 1; i++) {
            auto &c1 = generators[i];

            for (unsigned int j = i + 1; j < generators.size(); j++) {
                auto &c2 = generators[j];

                if (xq[i] * xq[j] >= 0)
                    continue;

                auto new_incidences = c1.incidences | c2.incidences;
                if (static_cast<int>(new_incidences.count()) > n + 1)
                    continue;

                bool is_redundant =
                    std::any_of(ret.generators.begin(), ret.generators.end(),
                                [&](Generator &c) {
                                    return c.incidences < new_incidences;
                                });
                if (is_redundant)
                    continue;

                rearranging_erase_if(ret.generators, [&](Generator &c) {
                    return new_incidences < c.incidences;
                });

                auto new_eq = intersect(c1.eq, c2.eq, x);
                if (new_eq.c > 0)
                    new_eq = -new_eq;

                ret.generators.push_back(Generator{new_eq, new_incidences});
            }
        }

        return ret;
    }

    bool average_if_nonzero(Hyperplane &q) const
    {
        q.c = 1;
        for (size_t i = 0; i < generators.size(); ++i)
            q.c = std::experimental::lcm(q.c, generators[i].eq.c);
        q.c = -q.c;

        for (size_t j = 0; j < dim; ++j) {
            q.a[j] = 0;
            for (size_t i = 0; i < generators.size(); ++i)
                q.a[j] += generators[i].eq.a[j] * (q.c / generators[i].eq.c);

            if (q.a[j] == 0)
                return false;
        }

        q.c *= static_cast<int>(generators.size());
        q.cancel();

        return true;
    }

    friend std::ostream &operator<<(std::ostream &os, const RationalCone &rhs)
    {
        for (auto &generator : rhs.generators)
            os << generator.eq << std::endl;
        return os;
    }
};

struct FinalCone {
    Hyperplane eq1, eq2;

    FinalCone() {}
    explicit FinalCone(const RationalCone &cone)
    {
        assert(cone.generators.size() == 2);

        std::array<std::array<Long, 2>, dim> blub1, blub2;

        for (int i = 0; i < dim; ++i) {
            blub1[i][0] = cone.generators[0].eq.a[i];
            blub1[i][1] = cone.generators[1].eq.a[i];

            blub2[i][1] = cone.generators[0].eq.a[i];
            blub2[i][0] = cone.generators[1].eq.a[i];
        }

        std::sort(blub1.begin(), blub1.end(),
                  [](auto &a, auto &b) { return a[0] > b[0]; });
        std::sort(blub2.begin(), blub2.end(),
                  [](auto &a, auto &b) { return a[0] > b[0]; });

        bool which = true;
        for (int i = 0; i < dim; ++i)
            if (blub1[i][0] != blub2[i][0])
                which = blub1[i][0] > blub2[i][0];

        auto &blub = which ? blub1 : blub2;

        eq1.c = 0;
        eq2.c = 0;
        for (int i = 0; i < dim; ++i) {
            eq1.a[i] = blub[i][0];
            eq2.a[i] = blub[i][1];
            eq1.c += eq1.a[i];
            eq2.c += eq2.a[i];
        }
        eq1.c = -eq1.c * 2 / two_times_r;
        eq2.c = -eq2.c * 2 / two_times_r;
    }

    operator RationalCone() const
    {
        RationalCone ret{};
        ret.generators.push_back({eq1, {}});
        ret.generators.push_back({eq2, {}});
        return ret;
    }

    bool operator<(const FinalCone &rhs) const
    {
        if (eq1 != rhs.eq1)
            return eq1 < rhs.eq1;

        if (eq2 != rhs.eq2)
            return eq2 < rhs.eq2;

        return false;
    }
};

std::set<FinalCone> deferred_cones;
int deferred_cones_insertions = 0;

std::set<RationalCone> sorted_q_cones;
int sorted_q_cones_insertions = 0;

class Stopwatch {
    std::chrono::time_point<std::chrono::high_resolution_clock> start;
public:
    Stopwatch() : start(std::chrono::high_resolution_clock::now()) {}

    double count() const {
        auto now = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = now - start;
        return duration.count();
    }
};

struct ClassificationData {
    template <class T, size_t N>
    using array = std::array<T, N>;

    array<Vector, dim + 1> x; // List of points that have to be
                              // allowed by the weight system
    array<array<Long, dim>, dim + 1> x_inner_q; // First index: n,
                                                // second index: number of qs
    array<Hyperplane, dim> qs;
    int f0[dim]; // TODO: This is something to check for redundant points

    array<RationalCone, dim> q_cones;

    weight_system_store_t *wli;
    Long wnum;    // Number of weight system candidates
    Long candnum; // Number of weight system candidates, including duplicates
    Long winum;   // Number of IP weight systems
    size_t recursion_level_counts[dim];
    size_t weight_counts[dim];
    Stopwatch stopwatch;
};

void print_stats(const ClassificationData &X)
{
    printf("%7.2f: ", X.stopwatch.count());
    for (size_t i = 1; i < dim - 1; ++i)
        printf("%ld ", X.recursion_level_counts[i]);
    printf("-- ");

    for (size_t i = 1; i < dim; ++i)
        printf("%ld ", X.weight_counts[i]);
    printf("-- ");

    printf("%ld -- %d ", X.candnum, weight_system_store_size(X.wli));
    printf("\n");
    fflush(stdout);
}

Long eval_eq(const Equation &e, const Vector &v)
{
    Long p = e.c;
    for (int i = 0; i < dim; ++i)
        p += v[i] * e.a[i];
    return p;
}

// Adds weight system wn to the sorted list X.wli, if it is basic and not
// already there.
void RgcAddweight(Hyperplane wn, ClassificationData &X)
{
    int i, j, p;

    if (!allow_weight_one_half)
        for (i = 0; i < dim; i++)
            if (2 * wn.a[i] == -wn.c)
                return;

    if (!allow_weight_one)
        for (i = 0; i < dim; i++)
            if (wn.a[i] == -wn.c)
                return;

    if (!allow_weights_sum_one)
        for (i = 0; i < dim - 1; i++)
            for (j = i + 1; j < dim; j++)
                if (wn.a[i] + wn.a[j] == -wn.c)
                    return;

    // Sort weights
    X.candnum++;
    for (i = 0; i < dim - 1; i++)
        for (p = i + 1; p < dim; p++)
            if (wn.a[i] > wn.a[p]) {
                Long k = wn.a[i];
                wn.a[i] = wn.a[p];
                wn.a[p] = k;
            } /* make n0<=n1<=...<=n# */

    // Add weight system to the sorted list if it is not already there
    weight_system_store_insert(X.wli, &wn);
}

int lex_cmp(const Vector &y1, const Vector &y2)
{
    for (int i = 0; i < dim; ++i) {
        if (y1[i] < y2[i])
            return -1;
        if (y1[i] > y2[i])
            return 1;
    }
    return 0;
}

// Returns true if the point should be trivially excluded
bool point_trivially_forbidden(const Vector &x, int ordered_from_coordinate)
{
    Long xsum = 0, xmax = 0;

    for (int l = 0; l < dim; l++) {
        xsum += x[l];
        if (x[l] > xmax)
            xmax = x[l];
    }

    if (!allow_weight_one && xsum == 1)
        return true;

    if (xsum == 2) {
        if (!allow_weight_one_half && xmax == 2)
            return true;

        if (!allow_weights_sum_one && xmax != 2)
            return true;
    }

    // Point does not allow positive weight systems (except if all coordinates
    // are 1)
    if (two_times_r == 2)
        if (xmax < 2)
            return true;

    // TODO: Why can we exclude this?
    if (two_times_r == 1)
        if (xmax < 3)
            return true;

    // TODO: This drastically removes redundant points. How does it work?
    if (!debug_ignore_symmetries)
        for (int l = ordered_from_coordinate; l < dim - 1; l++)
            if (x[l] < x[l + 1])
                return true;

    return false;
}

// Returns true if the point should be excluded
bool point_forbidden(const Vector &x, int n, ClassificationData &X)
{
    for (int i = 0; i < n - 1; ++i) {
        Long rel = X.x_inner_q[i + 1][i] - X.x_inner_q[n][i];
        if (rel > 0)
            return true;
        if (!disable_lex_cmp && rel == 0 && lex_cmp(X.x[i + 1], x) < 0)
            return true;
    }

    for (int i = 1; i < n; ++i) {
        auto &x_other = X.x[i];
        Long x_diff[dim];

        for (int j = 0; j < dim; ++j)
            x_diff[j] = x_other[j] - x[j];

        Long v = std::abs(x_diff[0]);

        for (int j = 1; j < dim; ++j) {
            if (x_diff[j] > 0)
                v = Fgcd(v, x_diff[j]);
            else if (x_diff[j] < 0)
                v = Fgcd(v, -x_diff[j]);
        }

        if (v != 1)
            for (int j = 0; j < dim; ++j)
                x_diff[j] /= v;

        bool all_positive = true;
        for (int j = 0; j < dim; ++j) {
            if (x_other[j] + x_diff[j] < 0) {
                all_positive = false;
                break;
            }
        }
        if (all_positive)
            return true;

        all_positive = true;
        for (int j = 0; j < dim; ++j) {
            if (x[j] - x_diff[j] < 0) {
                all_positive = false;
                break;
            }
        }
        if (all_positive)
            return true;
    }

    return false;
}

bool intersect_if_positive(Hyperplane &q, const Hyperplane &q1,
                           const Hyperplane &q2, const Vector &v)
{
    Long e1 = v * q1;
    if (e1 == 0)
        return false;

    Long e2 = v * q2;
    if (e1 < 0) {
        if (e2 <= 0)
            return false;
        e1 *= -1;
    } else {
        if (e2 >= 0)
            return false;
        e2 *= -1;
    }

    Long gcd = std::experimental::gcd(e1, e2);
    e1 /= gcd;
    e2 /= gcd;

    q.a = e2 * q1.a + e1 * q2.a;
    q.c = e2 * q1.c + e1 * q2.c;

    q.cancel();
    return true;
}

void RecConstructRgcWeights(int n, ClassificationData &X)
{
    /* we have q[n-1], x[n] */
    Hyperplane &q = X.qs[n];
    RationalCone &q_cone = X.q_cones[n];

    ++X.recursion_level_counts[n];

    switch (n) {
    case 0:
        X.f0[0] = 0;

        if (two_times_r % 2 == 0)
            q_cone = RationalCone::positive_orthant(two_times_r / 2, -1);
        else
            q_cone = RationalCone::positive_orthant(two_times_r, -2);

        if (!q_cone.average_if_nonzero(q))
            return;
        break;
    case dim - 1:
        assert(X.q_cones[dim - 2].generators.size() == 2);

        ++X.weight_counts[dim - 1];
        // if (++X.weight_counts[dim - 1] % 100000 == 0)
        //     print_stats(X);

        if (!intersect_if_positive(q, X.q_cones[dim - 2].generators[0].eq,
                                   X.q_cones[dim - 2].generators[1].eq,
                                   X.x[dim - 1]))
            return;
        break;
    default:
        auto &x = X.x[n];
        int i;
        for (i = dim - 1; (i >= X.f0[n - 1]) && (x[i] == 0); i--)
            ;
        X.f0[n] = ++i;

        q_cone = X.q_cones[n - 1].restrict(X.x[n], n);

        // sorted_q_cones.insert(X.q_cones[n]);
        // ++sorted_q_cones_insertions;

        // if (sorted_q_cones_insertions % 10000 == 0)
        //     printf("sorted_q_cones: %d/%d\n", (int)sorted_q_cones.size(),
        //     sorted_q_cones_insertions);

        if (!q_cone.average_if_nonzero(q))
            return;
    }

    // print_stats(x);

    if (n < dim - redundancy_check_skip_recursions) {
        if (n == dim - 1)
            q_cone = X.q_cones[n - 1].restrict(X.x[n], n);

        bool skip = false;

        enumerate_points_on(q, [&](auto &x) {
            for (unsigned int i = 1; i < q_cone.generators.size(); ++i) {
                if (q_cone.generators[i].eq * x != 0)
                    return;
            }

            // TODO: Do we have to check if x is linearly independent of the
            // other points?
            for (int i = 0; i < n; ++i) {
                Long diff = X.x_inner_q[i + 1][i] - x * X.qs[i];
                if (diff > 0 || (!disable_lex_cmp && diff == 0 && lex_cmp(x, X.x[i + 1]) > 0)
                    ) {
                    skip = true;
                    return;
                }
            }
        });

        if (skip)
            return;
    }

    RgcAddweight(q, X);

    if (defer_last_recursion && n == dim - 2) {
        ++deferred_cones_insertions;
        deferred_cones.insert(FinalCone{q_cone});
        if (deferred_cones_insertions % 100000 == 0) {
            printf("%7.2f: ", X.stopwatch.count());
            cout << deferred_cones.size() << "/" << deferred_cones_insertions
                 << std::endl;
        }
        return;
    }

    if (n == dim - 1)
        return;

    enumerate_points_below(q, [&](auto &x) {
        if (point_trivially_forbidden(x, X.f0[n]))
            return;

        for (int i = 0; i < n + 1; ++i)
            X.x_inner_q[n + 1][i] = X.qs[i] * x;

        if (point_forbidden(x, n + 1, X))
            return;

        X.x[n + 1] = x;
        RecConstructRgcWeights(n + 1, X);
    });
}

void AddPointToPoly(Long y[dim], PolyPointList &P)
{
    if (P.np >= POINT_Nmax) {
        printf("ohh no\n");
        exit(1);
    }

    int i;
    for (i = 0; i < P.n; i++)
        P.x[P.np][i] = y[i];
    P.np++;
}

int WsIpCheck(const Hyperplane &q)
{
    int k, l;
    auto P = std::unique_ptr<PolyPointList>(
        static_cast<PolyPointList *>(operator new(sizeof(PolyPointList))));
    VertexNumList V;
    EqList E;
    Vector y;
    Long yq[dim];
    P->n = dim;
    P->np = 0;
    for (k = 0; k < dim; k++) {
        y[k] = 0;
        yq[k] = 0;
    }
    k = dim - 1;
    AddPointToPoly(y.data(), *P);
    y[k] = -1; /* starting point just outside                       */
    yq[k] = -q.a[k];
    while (k >= 0) {
        y[k]++;
        yq[k] += q.a[k];
        for (l = k + 1; l < dim; l++)
            yq[l] = yq[k];
        if (yq[k] == -q.c)
            AddPointToPoly(y.data(), *P);
        for (k = dim - 1; (k >= 0 ? (yq[k] + q.a[k] > -q.c) : 0); k--)
            y[k] = 0;
    }
    /* sets k to the highest value where y[k] didn't exceed max value;
       resets the following max values to min values                 */
    if (P->np <= dim) {
        return 0;
    }
    Find_Equations(&*P, &V, &E);
    if (E.ne < dim) {
        return 0;
    }
    for (k = 0; k < dim; k++)
        y[k] = 1;
    for (k = 0; k < E.ne; k++)
        if (eval_eq(E.e[k], y) <= 0)
            if (!E.e[k].c) {
                return 0;
            }
    k = P->np - 1;
    return k;
}

void RgcWeights(void)
{
    ClassificationData X{};

    X.wli = weight_system_store_new();

    std::ofstream out{"out", std::ofstream::binary};
    std::ofstream txt_out{"out.txt"};

    if (read_from_file) {
        RationalCone cone{};
        cone.generators.resize(2);
        int in_count = 0;

        std::ifstream cones_in{"deferred_cones", std::ifstream::binary};
        srand(static_cast<unsigned>(time(NULL)));
        while (cones_in) {
            // cones_in.seekg((rand() % 46890549) * 20);
            Hyperplane q;

            q.c = 0;
            for (int i = 0; i < dim; ++i) {
                uint16_t v16;
                cones_in.read(reinterpret_cast<char *>(&v16), sizeof(v16));
                q.a[i] = ntohs(v16);
                q.c += q.a[i];
            }
            q.c = -q.c * 2 / two_times_r;
            cone.generators[0].eq = q;

            q.c = 0;
            for (int i = 0; i < dim; ++i) {
                uint16_t v16;
                cones_in.read(reinterpret_cast<char *>(&v16), sizeof(v16));
                q.a[i] = ntohs(v16);
                q.c += q.a[i];
            }
            q.c = -q.c * 2 / two_times_r;
            cone.generators[1].eq = q;

            if (!cone.average_if_nonzero(q)) {
                cout << ".";
                continue;
            }

            enumerate_points_below(q, [&](auto &x) {
                Hyperplane q_final;
                // if (point_trivially_forbidden(x, INT_MAX))
                //     return;

                if (intersect_if_positive(q_final, cone.generators[0].eq,
                                          cone.generators[1].eq, x)) {
                    // for (int i = 0; i < dim; ++i) {
                    //     auto v = q_final.a[i];
                    //     assert(v >= 0 && v <= UINT32_MAX);
                    //     uint32_t v32 = htonl(v);
                    //     out.write(reinterpret_cast<const char *>(&v32), sizeof(v32));
                    // }

                    RgcAddweight(q_final, X);
                }
            });

            ++in_count;

            // printf("%7.2f %d %d->%d\n", (time(NULL) - X.start_time) / 60.0,
            //        weight_system_store_size(X.wli), in_count, X.candnum);

            // if (weight_system_store_size(X.wli) >= 100000)
            //     break;
        }

        print_stats(X);
    } else {

    RecConstructRgcWeights(0, X);
    // fprintf(stderr, "\n");

    printf("q_cones: %ld/%d also %ld\n", sorted_q_cones.size(),
           sorted_q_cones_insertions, deferred_cones.size());

    print_stats(X);

    if (write_cones) {
        std::ofstream cones_out{"deferred_cones", std::ofstream::binary};
        for (auto &cone : deferred_cones) {
            for (int i = 0; i < dim; ++i) {
                auto v = cone.eq1.a[i];
                assert(v >= 0 && v <= UINT16_MAX);
                uint16_t v16 = htons(static_cast<uint16_t>(v));
                cones_out.write(reinterpret_cast<const char *>(&v16), sizeof(v16));
            }
            for (int i = 0; i < dim; ++i) {
                auto v = cone.eq2.a[i];
                assert(v >= 0 && v <= UINT16_MAX);
                uint16_t v16 = htons(static_cast<uint16_t>(v));
                cones_out.write(reinterpret_cast<const char *>(&v16), sizeof(v16));
            }
        }
        print_stats(X);
        return;
    }

    if (defer_last_recursion) {
        for (auto &cone : deferred_cones) {
            RationalCone q_cone = cone;
            Hyperplane q;
            if (!q_cone.average_if_nonzero(q))
                return;

            enumerate_points_below(q, [&](auto &x) {
                    Hyperplane q_final;
                    // if (point_trivially_forbidden(x, INT_MAX))
                    //     return;

                    if (intersect_if_positive(q_final, cone.eq1, cone.eq2, x))
                        RgcAddweight(q_final, X);
                });
        }

        print_stats(X);
    }

    }

    X.wnum = weight_system_store_size(X.wli);

    fflush(stdout);
    // for (auto &cone : sorted_q_cones)
    //     cout << cone << endl;

    // std::ofstream out{"out", std::ofstream::binary};
    // std::ofstream txt_out{"out.txt"};

    const Hyperplane *e;
    weight_system_store_begin_iteration(X.wli);
    while ((e = weight_system_store_next(X.wli))) {
        // for (int i = 0; i < dim; ++i)
        //     out.write(reinterpret_cast<const char *>(&e->a[i]), sizeof(e->a[i]));

        // txt_out << "{" << e->a[0];
        // for (int i = 1; i < dim; ++i)
        //     txt_out << "," << e->a[i];
        // txt_out << "},\n";

        int i = WsIpCheck(*e);
        if (i) {
            // cout << *e << endl;
            // printf("  np=%d\n", i);
            X.winum++;
            /*else PrintEquation(X.wli[i], dim, "n");*/
        }
    }
    print_stats(X);
    printf("#ip=%ld, #cand=%ld(%ld)\n", X.winum, X.wnum, X.candnum);
}
}

FILE *inFILE;
FILE *outFILE;

int main()
{
    inFILE = stdin;
    outFILE = stdout;

    RgcWeights();
    return 0;
}
