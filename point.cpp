#include "point.h"

bool is_sorted(const Point &x,
               const std::vector<std::pair<unsigned, unsigned>> &checks)
{
    for (const auto &check : checks)
        if (x.coords[check.first] < x.coords[check.second])
            return false;
    return true;
}
