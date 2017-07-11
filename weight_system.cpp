#include "weight_system.h"
#include <algorithm>
#include <experimental/numeric>

Long distance(const WeightSystem &ws, const Point &x)
{
    Long ret = 0;
    for (unsigned i = 0; i < dim; ++i)
        ret += (x.coords[i] * r_numerator - r_denominator) * ws.weights[i];
    return ret;
}

void cancel(WeightSystem &ws)
{
    if (dim == 0)
        return;

    Long gcd = std::abs(ws.weights[0]);

    for (unsigned i = 1; i < dim; ++i)
        gcd = std::experimental::gcd(gcd, ws.weights[i]);

    if (gcd != 1 && gcd != 0)
        ws /= gcd;
}

void sort(WeightSystem &ws)
{
    std::sort(ws.weights.begin(), ws.weights.end());
}

Long norm(const WeightSystem &ws)
{
    return std::accumulate(ws.weights.begin(), ws.weights.end(), 0);
}

WeightSystem intersect(const WeightSystem &q1, const WeightSystem &q2,
                       const Point &x)
{
    Long e1 = distance(q1, x);
    Long e2 = distance(q2, x);

    // Long gcd = std::experimental::gcd(e1, e2);
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
    int k = dim - 1;
    while (ax[k] + q.weights[k] * r_numerator >= 0) {
        if (k == 0)
            return false;
        x.coords[k] = 0;
        --k;
    }

    x.coords[k]++;
    ax[k] += q.weights[k] * r_numerator;
    for (int i = k + 1; i < dim; ++i)
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
        int k = dim - 1;
        while (ax[k] + q.weights[k] * r_numerator > 0) {
            if (k == 0)
                return false;
            x.coords[k] = 0;
            --k;
        }

        x.coords[k]++;
        ax[k] += q.weights[k] * r_numerator;
        for (int i = k + 1; i < dim; ++i)
            ax[i] = ax[k];

        if (ax[k] == 0)
            return true;
    }
}
