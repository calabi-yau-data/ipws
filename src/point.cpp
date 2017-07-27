#include "point.h"

using std::array;
using std::pair;
using std::vector;

namespace {

bool points_have_symmetry(gsl::span<const Point> points, unsigned idx1,
                          unsigned idx2)
{
    for (const auto &p : points)
        if (p.coords[idx1] != p.coords[idx2])
            return false;
    return true;
}
}

// the indices in the pairs are in ascending orders
std::vector<std::pair<unsigned, unsigned>> symmetries(
    gsl::span<const Point> points)
{
    vector<pair<unsigned, unsigned>> ret{};

    array<bool, dim> done{};
    done.fill(false);

    for (unsigned i = 0; i < dim - 1; ++i) {
        if (done[i])
            continue;

        unsigned k = i;
        for (unsigned j = i + 1; j < dim; ++j) {
            if (done[j])
                continue;

            if (points_have_symmetry(points, i, j)) {
                ret.push_back(pair<unsigned, unsigned>(k, j));
                k = j;
                done[j] = true;
            }
        }
    }

    return ret;
}

bool is_sorted(const Point &x,
               const std::vector<std::pair<unsigned, unsigned>> &checks)
{
    for (const auto &check : checks)
        if (x.coords[check.first] < x.coords[check.second])
            return false;
    return true;
}
