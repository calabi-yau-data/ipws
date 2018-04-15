#include "history.h"
#include "settings.h"

// TODO: verify this function
bool last_point_redundant(int n, const History &history)
{
    Point x = history.points[n];

    for (int i = 0; i < n; ++i) {
        Ring rel = history.point_weight_system_distances[n][i] -
                   history.point_weight_system_distances[i][i];
        if (rel < 0) {
            // std::cout << "x" << n + 1 << ".q" << i << " < x" << i + 1 << ".q"
            //           << i << "; x" << n + 1 << " = " << x << ", x" << i + 1
            //           << " = " << history.points[i] << ", q" << i << " = "
            //           << history.weight_systems[i] << std::endl;
            return true;
        }
        if (!g_settings.debug_disable_lex_order && rel == 0 &&
            x > history.points[i]) {
            // std::cout << "x" << n + 1 << ".q" << i << " = x" << i + 1 << ".q"
            //           << i << " and x" << n + 1 << " > x" << i + 1 << "; x"
            //           << n + 1 << " = " << x << ", x" << i + 1 << " = "
            //           << history.points[i] << ", q" << i << " = "
            //           << history.weight_systems[i] << std::endl;
            return true;
        }
    }

    for (int i = 0; i < n; ++i) {
        Point x_other = history.points[i];
        Point x_diff = x_other - x;

        Ring v = gcd(x_diff);
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
            Ring rel = distance(history.weight_systems[i], x) -
                       history.point_weight_system_distances[i][i];
            if (rel < 0) {
                // std::cout << "x.q" << i << " < x" << i + 1 << ".q" << i
                //           << "; x = " << x << ", x" << i + 1 << " = "
                //           << history.points[i] << ", q" << i << " = "
                //           << history.weight_systems[i] << std::endl;
                return true;
            }
            if (!g_settings.debug_disable_lex_order && rel == 0 &&
                x > history.points[i]) {
                // std::cout << "x.q" << i << " = x" << i + 1 << ".q" << i
                //           << " and x > x" << i + 1 << "; x = " << x << ", x"
                //           << i + 1 << " = " << history.points[i] << ", q" << i
                //           << " = " << history.weight_systems[i] << std::endl;
                return true;
            }
        }
    }

    return false;
}
