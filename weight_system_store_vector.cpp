#include "weight_system_store.h"

#include <algorithm>
#include <vector>

class WeightSystem {
    int16_t weights[DIMENSION];

public:
    bool operator<(const WeightSystem &rhs) const
    {
        for (size_t i = DIMENSION; i--;)
            if (weights[i] != rhs.weights[i])
                return weights[i] < rhs.weights[i];
        return false;
    }

    WeightSystem(const Equation &eq)
    {
        for (size_t i = 0; i < DIMENSION; ++i)
            weights[i] = eq.a[i];
    }

    operator Equation() const
    {
        Equation eq;

        eq.c = 0;
        for (size_t i = 0; i < DIMENSION; ++i) {
            eq.a[i] = weights[i];
            eq.c -= eq.a[i];
        }
        eq.c = eq.c * 2 / TWO_TIMES_R;

        return eq;
    }
};

using Collection = std::vector<WeightSystem>;

struct weight_system_store {
    Collection s;
    Collection::iterator it;
    Equation eq;
};

weight_system_store_t *weight_system_store_new()
{
    return new weight_system_store_t{};
}

void weight_system_store_free(weight_system_store_t *store)
{
    delete store;
}

void weight_system_store_insert(weight_system_store_t *store, const Equation *e)
{
    WeightSystem ws{*e};
    auto it = std::lower_bound(store->s.begin(), store->s.end(), ws);
    if (it == store->s.end() || ws < *it)
        store->s.insert(it, ws);
}

int weight_system_store_size(weight_system_store_t *store)
{
    return store->s.size();
}

void weight_system_store_begin_iteration(weight_system_store_t *store)
{
    store->it = store->s.begin();
}

const Equation *weight_system_store_next(weight_system_store_t *store)
{
    if (store->it == store->s.end())
        return nullptr;

    const WeightSystem &ws = *store->it++;
    store->eq = ws;

    return &store->eq;
}
