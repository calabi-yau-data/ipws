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

void checko(const WeightSystem &ws, Ring &picard_number,
            std::array<Ring, dim - 1> &hodge_numbers)
{
    using std::make_unique;

    assert(r_denominator == 1);
    assert(dim <= POLY_Dmax + 1);

    VertexNumList V;
    BaHo BH;

    // allocate large structures on the heap once
    static auto CW = make_unique<CWS>();
    static auto E = make_unique<EqList>();
    static auto DE = make_unique<EqList>();
    static auto P = make_unique<PolyPointList>();
    static auto DP = make_unique<PolyPointList>();
    static auto PM = make_unique<PairMatStruct>();
    static auto DPM = make_unique<PairMatStruct>();

    CW->N = dim;
    CW->index = r_numerator;
    CW->nw = 1;
    Ring norm = 0;
    for (ptrdiff_t i = 0; i < CW->N; i++) {
        CW->W[0][i] = ws.weights[i] * r_numerator;
        norm += ws.weights[i];
    }
    CW->d[0] = norm;
    CW->nz = 0;

    Make_CWS_Points(&*CW, &*P);

    bool reflexive = false;
    bool ip = Find_Equations(&*P, &V, &*E);
    if (!ip)
        return;

    reflexive = EL_to_PPL(&*E, &*DP, &P->n);

    // Sort_VL(&V); // unnecessary?

    Make_VEPM(&*P, &V, &*E, PM->pm);

    if (reflexive) {
        int ret = Transpose_PM(PM->pm, DPM->pm, V.nv, E->ne);
        assert(ret);
        VNL_to_DEL(&*P, &V, &*DE);
        Complete_Poly(DPM->pm, &*DE, E->ne, &*DP);
        RC_Calc_BaHo(&*P, &V, &*E, &*DP, &BH);
    } else {
        BH.mp = P->np;
        BH.mv = V.nv;
        BH.np = 0;
        BH.nv = E->ne;
    }

    assert(BH.n == dim - 1);

    for (size_t i = 0; i < hodge_numbers.size(); ++i)
        hodge_numbers[i] = BH.h1[i + 1];

    switch (dim) {
    case 4:
        picard_number = BH.pic;
        break;
    case 5:
        picard_number = 2 * (BH.h1[1] - BH.h1[2]);
        break;
    case 6:
        picard_number = 48 + 6 * (BH.h1[1] - BH.h1[2] + BH.h1[3]);
        break;
    default:
        assert(false);
    }
}
