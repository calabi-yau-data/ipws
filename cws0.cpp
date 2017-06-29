#include <time.h>
#include <algorithm>
#include <array>
#include <bitset>
#include <iostream>
#include <memory>
#include <set>
#include <vector>

#include "Global.h"
#include "Rat.h"

#include "weight_system_store.h"

namespace {

const bool allow11 = false;

struct EquationRedux {
    std::array<Long, dim> a;
    Long c;

    EquationRedux() {}
    EquationRedux(const Equation &e)
    {
        c = e.c;
        for (int i = 0; i < dim; ++i)
            a[i] = e.a[i];
    }

    bool operator!=(const EquationRedux &rhs) const
    {
        if (c != rhs.c)
            return true;

        for (int i = 0; i < dim; ++i)
            if (a[i] != rhs.a[i])
                return true;

        return false;
    }

    bool operator<(const EquationRedux &rhs) const
    {
        if (c != rhs.c)
            return c < rhs.c;

        for (int i = 0; i < dim; ++i)
            if (a[i] != rhs.a[i])
                return a[i] < rhs.a[i];

        return false;
    }

    friend std::ostream &operator<<(std::ostream &os, const EquationRedux &rhs)
    {
        os << rhs.c << "  ";

        for (auto x : rhs.a)
            os << x << " ";
        return os;
    }
};

class Cone {
    int ne;
    // std::array<EquationRedux, dim * 2> e;
    std::vector<EquationRedux> e;

public:
    Cone(const EqList &eq_list)
    {
        e.resize(eq_list.ne);

        if (eq_list.ne > (int)e.size()) {
            std::cout << "too much\n";
            exit(1);
        }

        ne = eq_list.ne;
        for (int i = 0; i < ne; ++i)
            e[i] = eq_list.e[i];

        std::sort(e.begin(), e.begin() + ne);
    }

    bool operator<(const Cone &rhs) const
    {
        if (ne != rhs.ne)
            return ne < rhs.ne;

        for (int i = 0; i < ne; ++i)
            if (e[i] != rhs.e[i])
                return e[i] < rhs.e[i];

        return false;
    }

    friend std::ostream &operator<<(std::ostream &os, const Cone &rhs)
    {
        for (int i = 0; i < rhs.ne; ++i)
            os << rhs.e[i] << std::endl;
        return os;
    }
};

// Enumerates all points x with nonnegative coordinates, satisfying
// q.a * x + q.c < 0.
template <class F>
void enumerate_points_below(const Equation &q, F callback)
{
    std::array<Long, dim> x;
    std::array<Long, dim> a_times_x;

    std::fill_n(x.begin(), dim - 1, 0);
    std::fill_n(a_times_x.begin(), dim - 1, 0);

    x[dim - 1] = -1;
    a_times_x[dim - 1] = -q.a[dim - 1];

    for (int k = dim - 1; k >= 0;) {
        x[k]++;
        a_times_x[k] += q.a[k];
        for (int i = k + 1; i < dim; i++)
            a_times_x[i] = a_times_x[k];

        callback(x);

        for (k = dim - 1; k >= 0 && a_times_x[k] + q.a[k] >= -q.c; k--)
            x[k] = 0;
    }
}

std::set<Cone> q_cones;
int q_cones_insertions = 0;

#define SKEW 0

void RgcWeights(void);

using Vector = std::array<Long, dim>;

template <size_t N>
bool operator<(const std::bitset<N> &lhs, const std::bitset<N> &rhs)
{
    return (lhs | rhs) == rhs;
}

struct RationalCone {
    using Incidence = std::bitset<dim>;

    EqList generators;
    std::array<Incidence, EQUA_Nmax> incidences;
};

struct ClassificationData {
    template <class T, size_t N>
    using array = std::array<T, N>;

    array<Vector, dim + 1> x; // List of points that have to be
                                        // allowed by the weight system
    array<array<Long, dim>, dim + 1> x_inner_q; // First index: n,
                                                // second index: number of q_tilde
    array<Equation, dim> q_tilde;
    int f0[dim]; // TODO: This is something to check for redundant points

    array<RationalCone, dim> muh;

    weight_system_store_t *wli;
    Long wnum;    // Number of weight system candidates
    Long candnum; // Number of weight system candidates, including duplicates
    Long winum;   // Number of IP weight systems
    size_t recursion_level_counts[dim];
    size_t weight_counts[dim];
    time_t start_time;
};

int ComputeAverageWeight(Equation &q, int n, ClassificationData &X);

void print_stats(const ClassificationData &X)
{
    printf("%8.3f: ", (time(NULL) - X.start_time) / 60.0);
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
void RgcAddweight(Equation wn, ClassificationData &X)
{
    int i, j, p, k;

    // Skip weight systems containing a weight of 1/2
    for (i = 0; i < dim; i++)
        if (2 * wn.a[i] == -wn.c)
            return;

    // Skip weight systems containing a weight of 1
    for (i = 0; i < dim; i++)
        if (wn.a[i] == -wn.c)
            return;

    // Skip weight systems containing two weights with a sum of 1
    if (!allow11)
        for (i = 0; i < dim - 1; i++)
            for (j = i + 1; j < dim; j++)
                if (wn.a[i] + wn.a[j] == -wn.c)
                    return;

    // Sort weights
    X.candnum++;
    for (i = 0; i < dim - 1; i++)
        for (p = i + 1; p < dim; p++)
            if (wn.a[i] > wn.a[p]) {
                k = wn.a[i];
                wn.a[i] = wn.a[p];
                wn.a[p] = k;
            } /* make n0<=n1<=...<=n# */

    // Add weight system to the sorted list if it is not already there
    weight_system_store_insert(X.wli, &wn);
}

void PrintQ(int n, ClassificationData &X)
{
    int i, j;
    assert(n < dim);
    for (i = 0; i < n; i++)
        printf(" ");
    printf("q: ne=%d\n", X.muh[n].generators.ne);
    for (j = 0; j < X.muh[n].generators.ne; j++) {
        printf("%d  ", (int)X.muh[n].generators.e[j].c);
        for (i = 0; i < dim; i++)
            printf("%d ", (int)X.muh[n].generators.e[j].a[i]);
        printf("\n");
    }
}

void PrintEquation(const Equation &q /*, char *c, int j*/)
{
    int i;
    printf("%d  ", (int)-q.c);
    for (i = 0; i < dim; i++)
        printf("%d ", (int)q.a[i]);
    /*printf("  %s  np=%d\n", c, j);*/
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

// Tests if the point y should be considered
int PointForbidden(const Vector &y, int n, ClassificationData &X)
{
    int l;
    Long ysum = 0, ymax = 0;

    // assert(n < dim);

    for (l = 0; l < dim; l++) {
        ysum += y[l];
        if (y[l] > ymax)
            ymax = y[l];
    }

    // Point leads to weight systems containing a weight of 1
    if (ysum < 2)
        return 1;

    // Point leads to weight systems containing a weight of 1/2 or two weights
    // with a sum of 1
    if (ysum == 2)
        if (!allow11 || ymax == 2)
            return 1;

    // Point does not allow positive weight systems (except if all coordinates
    // are 1)
    if (two_times_r == 2)
        if (ymax < 2)
            return 1;

    // TODO: Why can we exclude this?
    if (two_times_r == 1)
        if (ymax < 3)
            return 1;

#if SKEW == 0
    // TODO: This drastically removes redundant points. How does it work?
    for (l = X.f0[n - 1]; l < dim - 1; l++)
        if (y[l] < y[l + 1])
            return 1;
#endif

    for (int i = 0; i < n; ++i)
        X.x_inner_q[n][i] = eval_eq(X.q_tilde[i], y);

    for (int i = 0; i < n - 1; ++i) {
        int rel = X.x_inner_q[i + 1][i] - X.x_inner_q[n][i];
        if (rel > 0)
            return 1;
        if (rel == 0 && lex_cmp(X.x[i + 1], y) < 0)
            return 1;
    }

    for (int i = 1; i < n; ++i) {
        auto &y_other = X.x[i];
        Long y_diff[dim];

        for (int j = 0; j < dim; ++j)
            y_diff[j] = y_other[j] - y[j];

        int v = 1;
        if (y_diff[0] > 0)
            v = y_diff[0];
        else if (y_diff[0] < 0)
            v = -y_diff[0];

        for (int j = 1; j < dim; ++j) {
            if (y_diff[j] > 0)
                v = Fgcd(v, y_diff[j]);
            else if (y_diff[j] < 0)
                v = Fgcd(v, -y_diff[j]);
        }

        if (v != 1)
            for (int j = 0; j < dim; ++j)
                y_diff[j] /= v;

        bool all_positive = true;
        for (int j = 0; j < dim; ++j) {
            if (y_other[j] + y_diff[j] < 0) {
                all_positive = false;
                break;
            }
        }
        if (all_positive)
            return 1;

        all_positive = true;
        for (int j = 0; j < dim; ++j) {
            if (y[j] - y_diff[j] < 0) {
                all_positive = false;
                break;
            }
        }
        if (all_positive)
            return 1;
    }

    return 0;
}

// Initializes X.q such that it represents the following weight systems:
// (r, 0, ... 0), (0, r, ... 0), ... (0, 0, ... r).
// TODO: Also does something to X.f0 and X.qI
bool ComputeQ0(Equation &q, ClassificationData &X)
{
    X.muh[0].generators.ne = dim;
    for (int i = 0; i < X.muh[0].generators.ne; i++) {
        for (int j = 0; j < dim; j++)
            X.muh[0].generators.e[i].a[j] = 0;
        if (two_times_r % 2) {
            X.muh[0].generators.e[i].a[i] = two_times_r;
            X.muh[0].generators.e[i].c = -2;
        } else {
            X.muh[0].generators.e[i].a[i] = two_times_r / 2;
            X.muh[0].generators.e[i].c = -1;
        }
    }
    X.f0[0] = 0;
    for (int i = 0; i < X.muh[0].generators.ne; i++)
        X.muh[0].incidences[i].set(i);

    return ComputeAverageWeight(q, 0, X);
}

bool IsRedundant(const RationalCone::Incidence &new_incidence,
                 const std::array<RationalCone::Incidence, EQUA_Nmax> &qINew,
                 int ne)
{
    for (int i = 0; i < ne; i++)
        if (qINew[i] < new_incidence)
            return true;
    return false;
}

bool ComputeQ(Equation &q, int n, ClassificationData &X)
{
    /* q[n] from q[n-1], x[n] */
    int i, j, k;
    auto &y = X.x[n];
    Long yqOld[EQUA_Nmax];
    EqList &qOld = X.muh[n - 1].generators;
    EqList &qNew = X.muh[n].generators;
    auto &qIOld = X.muh[n - 1].incidences;
    auto &qINew = X.muh[n].incidences;

    assert(n < dim);
    for (i = dim - 1; (i >= X.f0[n - 1]) && (y[i] == 0); i--)
        ;
    X.f0[n] = ++i;

    qNew.ne = 0;
    for (i = 0; i < qOld.ne; i++)
        if (!(yqOld[i] = eval_eq(qOld.e[i], y))) {
            qNew.e[qNew.ne] = qOld.e[i];
            qINew[qNew.ne] = qIOld[i];
            (qNew.ne)++;
        }

    for (i = 0; i < qOld.ne - 1; i++)
        for (j = i + 1; j < qOld.ne; j++)
            if (yqOld[i] * yqOld[j] < 0) {
                auto new_incidence = qIOld[i] | qIOld[j];

                if (static_cast<int>(new_incidence.count()) > n + 1)
                    continue;

                if (IsRedundant(new_incidence, qINew, qNew.ne))
                    continue;

                for (k = qNew.ne - 1; k >= 0; k--)
                    if (new_incidence < qINew[k]) {
                        qINew[k] = qINew[qNew.ne - 1];
                        qNew.e[k] = qNew.e[qNew.ne - 1];
                        qNew.ne--;
                    }
                assert(qNew.ne < EQUA_Nmax - 1);
                qINew[qNew.ne] = new_incidence;
                qNew.e[qNew.ne] =
                    EEV_To_Equation(&qOld.e[i], &qOld.e[j], y.data(), dim);
                if (qNew.e[qNew.ne].c > 0) {
                    for (k = 0; k < dim; k++)
                        qNew.e[qNew.ne].a[k] *= -1;
                    qNew.e[qNew.ne].c *= -1;
                }
                qNew.ne++;
            }

    q_cones.insert(qNew);
    ++q_cones_insertions;

    if (q_cones_insertions % 10000 == 0)
        printf("q_cones: %d/%d\n", (int)q_cones.size(), q_cones_insertions);

    return ComputeAverageWeight(q, n, X);
}

Long Flcm(Long a, Long b)
{
    return (a * b) / Fgcd(a, b);
}

// Normalizes q such that it does not have any common divisors
void Cancel(Equation &q)
{
    Long gcd = -q.c;
    int j;
    for (j = 0; j < dim; j++)
        gcd = Fgcd(gcd, q.a[j]);
    if (gcd > 1) {
        for (j = 0; j < dim; j++)
            q.a[j] /= gcd;
        q.c /= gcd;
    }
}

int ComputeAverageWeight(Equation &q, int n, ClassificationData &X)
{
    int i, j;
    EqList &el = X.muh[n].generators;

    ++X.weight_counts[n];

    if (el.ne < dim - n)
        return 0;

    q.c = -1;
    for (i = 0; i < el.ne; i++) {
        if (el.e[i].c >= 0) {
            PrintQ(n, X);
            exit(0);
        }
        q.c = -Flcm(-q.c, -el.e[i].c);
    }
    for (j = 0; j < dim; j++) {
        q.a[j] = 0;
        for (i = 0; i < el.ne; i++)
            q.a[j] += el.e[i].a[j] * (q.c / el.e[i].c);

        q.a[j] += el.e[0].a[j] * (q.c / el.e[0].c) * SKEW;

        if (q.a[j] <= 0) {
            assert(q.a[j] == 0);
            return 0;
        }
    }

    q.c *= el.ne + SKEW;
    Cancel(q);
    return 1;
}

bool ComputeLastQ(Equation &q, ClassificationData &X)
{
    /* q[dim-1] from q[dim-2], x[dim-1] */
    assert(X.muh[dim - 2].generators.ne == 2);

    int i;
    Equation &q0 = X.muh[dim - 2].generators.e[0];
    Equation &q1 = X.muh[dim - 2].generators.e[1];
    auto &y = X.x[dim - 1];

    ++X.weight_counts[dim - 1];
    // if (++X.weight_counts[dim - 1] % 100000 == 0)
    //     print_stats(X);

    Long yq0 = eval_eq(q0, y);
    if (!yq0)
        return false;

    Long yq1 = eval_eq(q1, y);
    if (yq0 < 0) {
        if (yq1 <= 0)
            return false;
        yq0 *= -1;
    } else {
        if (yq1 >= 0)
            return false;
        yq1 *= -1;
    }
    q.c = yq0 * q1.c + yq1 * q0.c;
    for (i = 0; i < dim; i++)
        q.a[i] = yq0 * q1.a[i] + yq1 * q0.a[i];
    Cancel(q);
    return true;
}

void RecConstructRgcWeights(int n, ClassificationData &X)
{
    /* we have q[n-1], x[n] */
    Equation q;

    ++X.recursion_level_counts[n];

    switch (n) {
    case 0:
        if (!ComputeQ0(q, X))
            return;
        break;
    case dim - 1:
        if (!ComputeLastQ(q, X))
            return;
        break;
    default:
        if (!ComputeQ(q, n, X))
            return;
    }

    // print_stats(x);

    RgcAddweight(q, X);

    if (n == dim - 1)
        return;

    X.q_tilde[n] = q;

    enumerate_points_below(q, [&](auto &x) {
        if (PointForbidden(x, n + 1, X))
            return;

        std::copy(x.begin(), x.end(), std::begin(X.x[n + 1]));
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

int WsIpCheck(const Equation &q)
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
    X.start_time = time(NULL);

    RecConstructRgcWeights(0, X);
    // fprintf(stderr, "\n");

    X.wnum = weight_system_store_size(X.wli);

    printf("q_cones: %d/%d\n", (int)q_cones.size(), q_cones_insertions);
    fflush(stdout);
    // for (auto &cone : q_cones)
    //     std::cout << cone << std::endl;

    print_stats(X);

    const Equation *e;
    weight_system_store_begin_iteration(X.wli);
    while (e = weight_system_store_next(X.wli)) {
        int i = WsIpCheck(*e);
        if (i) {
            // PrintEquation(*e);
            // printf("  np=%d\n", i);
            X.winum++;
            /*else PrintEquation(X.wli[i], dim, "n");*/
        }
    }
    printf("#ip=%ld, #cand=%ld(%ld)\n", X.winum, X.wnum, X.candnum);
}
};

FILE *inFILE;
FILE *outFILE;

int main()
{
    inFILE = stdin;
    outFILE = stdout;

    RgcWeights();
    return 0;
}
