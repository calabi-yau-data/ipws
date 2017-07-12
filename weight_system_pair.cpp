#include "weight_system_pair.h"

#include <array>

using std::array;
using std::pair;
using std::vector;

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
        ret.weights[i] =
            pair.first.weights[i] * norm2 + pair.second.weights[i] * norm1;

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

// namespace {
// bool has_symmetry(const WeightSystemPair &pair, unsigned idx1,
//                   unsigned idx2, bool permuted)
// {

// }
// }

// std::vector<std::pair<unsigned, unsigned>> symmetries(
//     const WeightSystemPair &pair)
// {
//     vector<std::pair<unsigned, unsigned>> ret{};
//     return ret;

//     array<bool, dim> done{};
//     done.fill(false);

//     for (unsigned i = 0; i < dim - 1; ++i) {
//         if (done[i])
//             continue;

//         unsigned k = i;
//         for (unsigned j = i + 1; j < dim; ++j) {
//             if (done[j])
//                 continue;

//             if (has_symmetry(pair, i, j, false) ||
//                 has_symmetry(pair, i, j, true)) {
//                 ret.push_back(std::pair<unsigned, unsigned>(k, j));
//                 k = j;
//                 done[j] = true;
//             }
//         }
//     }

//     return ret;
// }
