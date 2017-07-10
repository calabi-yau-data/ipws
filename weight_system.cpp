#include <algorithm>
#include <experimental/numeric>

#include "weight_system.h"

Long distance(const WeightSystem &ws, const Vector &x)
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
                       const Vector &x)
{
    Long e1 = distance(q1, x);
    Long e2 = distance(q2, x);

    // Long gcd = std::experimental::gcd(e1, e2);
    // e1 /= gcd;
    // e2 /= gcd;

    WeightSystem ret{};
    if (e1 < 0)
        ret = e1 * q2 - e2 * q1;
    else
        ret = e2 * q1 - e1 * q2;

    cancel(ret);
    return ret;
}
