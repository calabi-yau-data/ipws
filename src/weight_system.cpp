#include "weight_system.h"
#include <algorithm>
#include <boost/math/common_factor.hpp>
#include <gsl/gsl>
#include "palp/Global.h"
#include "settings.h"

using gsl::span;
using std::array;

Ring distance(const WeightSystem<dim> &ws, const Point &x)
{
    Ring ret = 0;
    for (unsigned i = 0; i < dim; ++i)
        ret += (x.coords[i] * r_numerator - r_denominator) * ws.weights[i];
    return ret;
}

void cancel(WeightSystem<dim> &ws)
{
    if (dim == 0)
        return;

    Ring gcd = std::abs(ws.weights[0]);

    for (unsigned i = 1; i < dim; ++i)
        gcd = boost::math::gcd(gcd, ws.weights[i]);

    if (gcd != 1 && gcd != 0)
        ws /= gcd;
}

void sort(WeightSystem<dim> &ws)
{
    std::sort(ws.weights.begin(), ws.weights.end());
}

Ring norm(const WeightSystem<dim> &ws)
{
    return std::accumulate(ws.weights.begin(), ws.weights.end(), 0);
}

const WeightSystem<dim> intersect(const WeightSystem<dim> &q1,
                                  const WeightSystem<dim> &q2, const Point &x)
{
    Ring e1 = distance(q1, x);
    Ring e2 = distance(q2, x);

    // Ring gcd = boost::math::gcd(e1, e2);
    // e1 /= gcd;
    // e2 /= gcd;

    WeightSystem<dim> ret{};
    if (e1 < 0)
        ret = e2 * q1 - e1 * q2;
    else
        ret = e1 * q2 - e2 * q1;

    cancel(ret);
    return ret;
}

WeightSystemPointsBelow::WeightSystemPointsBelow(const WeightSystem<dim> &q)
        : q{q}
{
    x.coords.fill(0);
    ax.fill(-norm(q) * r_denominator);

    x.coords[dim - 1] -= 1;
    ax[dim - 1] -= q.weights[dim - 1] * r_numerator;
}

bool WeightSystemPointsBelow::find_next()
{
    unsigned k = dim - 1;
    while (ax[k] + q.weights[k] * r_numerator >= 0) {
        if (k == 0)
            return false;
        x.coords[k] = 0;
        --k;
    }

    x.coords[k]++;
    ax[k] += q.weights[k] * r_numerator;
    for (unsigned i = k + 1; i < dim; ++i)
        ax[i] = ax[k];

    return true;
}

WeightSystemPointsOn::WeightSystemPointsOn(const WeightSystem<dim> &q) : q{q}
{
    x.coords.fill(0);
    ax.fill(-norm(q) * r_denominator);

    x.coords[dim - 1] -= 1;
    ax[dim - 1] -= q.weights[dim - 1] * r_numerator;
}

bool WeightSystemPointsOn::find_next()
{
    while (true) {
        unsigned k = dim - 1;
        while (ax[k] + q.weights[k] * r_numerator > 0) {
            if (k == 0)
                return false;
            x.coords[k] = 0;
            --k;
        }

        x.coords[k]++;
        ax[k] += q.weights[k] * r_numerator;
        for (unsigned i = k + 1; i < dim; ++i)
            ax[i] = ax[k];

        if (ax[k] == 0)
            return true;
    }
}

static void add_point(const Point &x, PolyPointList *P)
{
    assert(P->np < POINT_Nmax);
    for (unsigned i = 0; i < dim; i++)
        P->x[P->np][i] = x.coords[i];
    P->np++;
}

static PolyPointList *new_point_list(const WeightSystem<dim> &ws)
{
    PolyPointList *P = (PolyPointList *)malloc(sizeof(PolyPointList));
    assert(P != nullptr);

    P->np = 0;
    P->n = dim;

    add_point(Point{}, P);

    auto gen = WeightSystemPointsOn(ws);
    while (gen.find_next())
        add_point(gen.get(), P);

    return P;
}

bool has_ip(const WeightSystem<dim> &ws)
{
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

bool good_weight_system(const WeightSystem<dim> &ws)
{
    Ring n = norm(ws);

    for (unsigned i = 0; i < dim; ++i) {
        if (!g_settings.allow_weights_one_half &&
            2 * ws.weights[i] * r_numerator == n * r_denominator)
            return false;

        if (!g_settings.allow_weights_one &&
            ws.weights[i] * r_numerator == n * r_denominator)
            return false;
    }

    if (!g_settings.allow_weights_sum_one)
        for (unsigned i = 0; i < dim - 1; ++i)
            for (unsigned j = i + 1; j < dim; ++j)
                if ((ws.weights[i] + ws.weights[j]) * r_numerator ==
                    n * r_denominator)
                    return false;

    return true;
}

struct PairMatStruct {
    PairMat pm;
};

#define MAXLD (26)

typedef struct {
    long n_nonIP, n_IP_nonRef, n_ref, // numbers of WS of certain types
        max_w, nr_max_w, // maximum weight in the reflexive/non-reflexive cases
        nr_n_w[MAXLD], n_w[MAXLD]; // numbers of weights of given [ld]
    int nr_max_mp, nr_max_mv, nr_max_nv, max_mp, max_mv, max_np, max_nv,
        max_h22, max_h1[POLY_Dmax - 1], // max values of certain entries of BH
        min_chi, max_chi,
        max_nf[POLY_Dmax + 1]; // range for chi, max facet numbers
} C5stats;

int int_ld(Long w)
{
    int i = -1;
    while (w) {
        w /= 2;
        i++;
    }
    return i;
}

void Initialize_C5S(C5stats *_C5S, int n)
{
    int k;
    _C5S->n_nonIP = 0;
    _C5S->n_IP_nonRef = 0;
    _C5S->n_ref = 0;
    _C5S->nr_max_mp = 0;
    _C5S->nr_max_mv = 0;
    _C5S->nr_max_nv = 0;
    _C5S->nr_max_w = 0;
    for (k = 0; k < MAXLD; k++) {
        _C5S->n_w[k] = 0;
        _C5S->nr_n_w[k] = 0;
    }
    _C5S->max_mp = 0;
    _C5S->max_mv = 0;
    _C5S->max_np = 0;
    _C5S->max_nv = 0;
    _C5S->max_h22 = 0;
    _C5S->max_w = 0;
    for (k = 1; k < n - 1; k++)
        _C5S->max_h1[k] = 0;
    for (k = 0; k <= n; k++)
        _C5S->max_nf[k] = 0;
    _C5S->max_chi = -100000000;
    _C5S->min_chi = 100000000;
}

extern "C" int QuickAnalysis(PolyPointList *_P, BaHo *_BH, FaceInfo *_FI);

void checko(const WeightSystem<dim> &ws, Ring &picard_number,
            std::array<Ring, dim - 1> &hodge_numbers, bool &reflexive)
{
    using std::make_unique;

    assert(r_denominator == 1);
    assert(dim <= POLY_Dmax + 1);

    // allocate large structures on the heap once
    static auto CW = make_unique<CWS>();
    static auto P = make_unique<PolyPointList>();

    CW->N = dim;
    CW->index = r_numerator;
    CW->nw = 1;
    Ring norm = 0;
    for (size_t i = 0; i < dim; i++) {
        CW->W[0][i] = ws.weights[i] * r_numerator;
        norm += ws.weights[i];
    }
    CW->d[0] = norm;
    CW->nz = 0;

    Make_CWS_Points(&*CW, &*P);

    BaHo BH;
    FaceInfo FI;
    C5stats C5S;

    if (!QuickAnalysis(&*P, &BH, &FI)) { // non-IP
        C5S.n_nonIP++;
        return;
    }
    Print_CWH(&*CW, &BH);

    reflexive = BH.np > 0;

    if (reflexive) {
        int i, chi = 48 + 6 * (BH.h1[1] - BH.h1[2] + BH.h1[3]),
            ld = int_ld(CW->W[0][5]);
        assert(0 <= ld);
        assert(ld < MAXLD);
        if (BH.mp > C5S.max_mp)
            C5S.max_mp = BH.mp;
        if (BH.mv > C5S.max_mv)
            C5S.max_mv = BH.mv;
        if (BH.np > C5S.max_np)
            C5S.max_np = BH.np;
        if (BH.nv > C5S.max_nv)
            C5S.max_nv = BH.nv;
        if (BH.h22 > C5S.max_h22)
            C5S.max_h22 = BH.h22;
        for (i = 1; i < BH.n - 1; i++)
            if (BH.h1[i] > C5S.max_h1[i])
                C5S.max_h1[i] = BH.h1[i];
        for (i = 0; i <= BH.n; i++)
            if (FI.nf[i] > C5S.max_nf[i])
                C5S.max_nf[i] = FI.nf[i];
        if (chi > C5S.max_chi)
            C5S.max_chi = chi;
        if (chi < C5S.min_chi)
            C5S.min_chi = chi;
        C5S.n_ref++;
        C5S.n_w[ld]++;
        if (CW->W[0][5] > C5S.max_w)
            C5S.max_w = CW->W[0][5];
    } else {
        if (BH.mp > C5S.nr_max_mp)
            C5S.nr_max_mp = BH.mp;
        if (BH.mv > C5S.nr_max_mv)
            C5S.nr_max_mv = BH.mv;
        if (BH.nv > C5S.nr_max_nv)
            C5S.nr_max_nv = BH.nv;
        C5S.n_IP_nonRef++;
        C5S.nr_n_w[int_ld(CW->W[0][5])]++;
        if (CW->W[0][5] > C5S.nr_max_w)
            C5S.nr_max_w = CW->W[0][5];
    }

    printf("non-IP: #=%ld\n", C5S.n_nonIP);
    printf(
        "IP, non-reflexive: #=%ld, max_mp=%d, max_mv=%d, max_nv=%d, "
        "max_w=%ld\n",
        C5S.n_IP_nonRef, C5S.nr_max_mp, C5S.nr_max_mv, C5S.nr_max_nv,
        C5S.nr_max_w);
    printf("  #(w5) of given ld: ");
    for (int i = 0; i < MAXLD; i++)
        printf(" %d:%ld", i, C5S.nr_n_w[i]);
    puts("");
    printf(
        "reflexive: #=%ld, max_mp=%d, max_mv=%d, max_np=%d, max_nv=%d, "
        "max_w=%ld\n",
        C5S.n_ref, C5S.max_mp, C5S.max_mv, C5S.max_np, C5S.max_nv,
        C5S.max_w);
    printf("  #(w5) of given ld: ");
    for (int i = 0; i < MAXLD; i++)
        printf(" %d:%ld", i, C5S.n_w[i]);
    puts("");
    printf("  max #(faces): %d %d %d %d %d\n", C5S.max_nf[0], C5S.max_nf[1],
           C5S.max_nf[2], C5S.max_nf[3], C5S.max_nf[4]);
    printf("  h11<=%d, h12<=%d, h13<=%d, h22<=%d, %d<=chi<=%d\n",
           C5S.max_h1[1], C5S.max_h1[2], C5S.max_h1[3], C5S.max_h22,
           C5S.min_chi, C5S.max_chi);

    // assert(BH.n == dim - 1);

    // for (size_t i = 0; i < hodge_numbers.size(); ++i)
    //     hodge_numbers[i] = BH.h1[i + 1];

    // switch (dim) {
    // case 4:
    //     picard_number = BH.pic;
    //     break;
    // case 5:
    //     picard_number = 2 * (BH.h1[1] - BH.h1[2]);
    //     break;
    // case 6:
    //     picard_number = 48 + 6 * (BH.h1[1] - BH.h1[2] + BH.h1[3]);
    //     break;
    // default:
    //     assert(false);
    // }
}
