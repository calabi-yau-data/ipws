#include "weight_system_pair.h"

WeightSystemPair canonicalize(const WeightSystemPair &pair)
{
    WeightSystemPair ret{};

    std::array<std::array<Long, 2>, dim> option0, option1;

    for (unsigned i = 0; i < dim; ++i) {
        option0[i][0] = pair.first.weights[i];
        option0[i][1] = pair.second.weights[i];

        option1[i][1] = pair.first.weights[i];
        option1[i][0] = pair.second.weights[i];
    }

    auto cmp = [](auto &a, auto &b) {
        if (a[0] != b[0])
            return a[0] > b[0];
        return a[1] > b[1];
    };

    std::sort(option0.begin(), option0.end(), cmp);
    std::sort(option1.begin(), option1.end(), cmp);

    bool choose_option0 = false;
    for (unsigned i = 0; i < dim; ++i) {
        if (cmp(option0[i], option1[i])) {
            choose_option0 = true;
            break;
        }
    }

    auto &option = choose_option0 ? option0 : option1;

    for (unsigned i = 0; i < dim; ++i) {
        ret.first.weights[i] = option[i][0];
        ret.second.weights[i] = option[i][1];
    }

    return ret;
}

WeightSystem average(const WeightSystemPair &pair)
{
    Long norm1 = norm(pair.first);
    Long norm2 = norm(pair.second);

    WeightSystem ret{};

    for (unsigned i = 0; i < dim; ++i)
        ret.weights[i] = pair.first.weights[i] * norm2
            + pair.second.weights[i] * norm1;

    cancel(ret);
    return ret;
}

bool restrict(const WeightSystemPair &pair, const Point &x, WeightSystem &ws)
{
    Long e1 = distance(pair.first, x);
    Long e2 = distance(pair.second, x);

    if (e1 < 0) {
        if (e2 <= 0)
            return false;
        ws = e2 * pair.first - e1 * pair.second;
    } else if (e1 > 0) {
        if (e2 >= 0)
            return false;
        ws = e1 * pair.second - e2 * pair.first;
    } else {
        return false;
    }

    cancel(ws);

    return true;
}
