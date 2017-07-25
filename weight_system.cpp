#include "weight_system.h"
#include <algorithm>
#include <boost/math/common_factor.hpp>
#include <gsl/gsl>
#include "io.h"
#include "settings.h"

extern "C" {
#include "palp/Global.h"
}

using gsl::span;
using std::array;

Ring distance(const WeightSystem &ws, const Point &x)
{
    Ring ret = 0;
    for (unsigned i = 0; i < dim; ++i)
        ret += (x.coords[i] * r_numerator - r_denominator) * ws.weights[i];
    return ret;
}

void cancel(WeightSystem &ws)
{
    if (dim == 0)
        return;

    Ring gcd = std::abs(ws.weights[0]);

    for (unsigned i = 1; i < dim; ++i)
        gcd = boost::math::gcd(gcd, ws.weights[i]);

    if (gcd != 1 && gcd != 0)
        ws /= gcd;
}

void sort(WeightSystem &ws)
{
    std::sort(ws.weights.begin(), ws.weights.end());
}

Ring norm(const WeightSystem &ws)
{
    return std::accumulate(ws.weights.begin(), ws.weights.end(), 0);
}

const WeightSystem intersect(const WeightSystem &q1, const WeightSystem &q2,
                             const Point &x)
{
    Ring e1 = distance(q1, x);
    Ring e2 = distance(q2, x);

    // Ring gcd = boost::math::gcd(e1, e2);
    // e1 /= gcd;
    // e2 /= gcd;

    WeightSystem ret{};
    if (e1 < 0)
        ret = e2 * q1 - e1 * q2;
    else
        ret = e1 * q2 - e2 * q1;

    cancel(ret);
    return ret;
}

WeightSystemPointsBelow::WeightSystemPointsBelow(const WeightSystem &q) : q{q}
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

WeightSystemPointsOn::WeightSystemPointsOn(const WeightSystem &q) : q{q}
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

static PolyPointList *new_point_list(const WeightSystem &ws)
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

bool has_ip(const WeightSystem &ws)
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

bool good_weight_system(const WeightSystem &ws)
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

void read(std::istream &f, WeightSystem &ws)
{
    for (unsigned i = 0; i < dim; ++i) {
        auto v = read_varint(f);
        assert(v <= std::numeric_limits<Ring>::max());
        ws.weights[i] = static_cast<Ring>(v);
    }
}

void write(std::ostream &f, const WeightSystem &ws)
{
    for (unsigned i = 0; i < dim; ++i)
        write_varint(f, ws.weights[i]);
}
