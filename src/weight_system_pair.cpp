#include "weight_system_pair.h"

#include <array>

using std::array;
using std::pair;
using std::vector;

const WeightSystemPair canonicalize(const WeightSystemPair &pair)
{
    WeightSystemPair ret{};

    std::array<std::array<Ring, 2>, dim> option0, option1;

    for (unsigned i = 0; i < dim; ++i) {
        option0[i][0] = pair.first.weights[i];
        option0[i][1] = pair.second.weights[i];

        option1[i][1] = pair.first.weights[i];
        option1[i][0] = pair.second.weights[i];
    }

    std::sort(option0.begin(), option0.end());
    std::sort(option1.begin(), option1.end());

    auto &option = option0 > option1 ? option0 : option1;

    for (unsigned i = 0; i < dim; ++i) {
        ret.first.weights[i] = option[i][0];
        ret.second.weights[i] = option[i][1];
    }

    return ret;
}

const WeightSystem<dim> average(const WeightSystemPair &pair)
{
    Ring norm1 = norm(pair.first);
    Ring norm2 = norm(pair.second);

    WeightSystem<dim> ret{};

    for (unsigned i = 0; i < dim; ++i)
        ret.weights[i] =
            pair.first.weights[i] * norm2 + pair.second.weights[i] * norm1;

    cancel(ret);
    return ret;
}

bool restrict(const WeightSystemPair &pair, const Point &x,
              WeightSystem<dim> &ws)
{
    Ring e1 = distance(pair.first, x);
    Ring e2 = distance(pair.second, x);

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
