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
