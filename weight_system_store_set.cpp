#include "weight_system_store.h"

#include <set>

struct WeightSystem {
    int16_t weights[DIMENSION];

    bool operator<(const WeightSystem &rhs) const
    {
        for (size_t i = DIMENSION; i--;)
            if (weights[i] != rhs.weights[i])
                return weights[i] < rhs.weights[i];
        return false;
    }
};

using Set = std::set<WeightSystem>;

struct weight_system_store {
    Set s;
    Set::iterator it;
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
    WeightSystem ws;
    for (size_t i = 0; i < DIMENSION; ++i)
        ws.weights[i] = e->a[i];

    store->s.insert(ws);
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

    store->eq.c = 0;
    for (size_t i = 0; i < DIMENSION; ++i) {
        store->eq.a[i] = ws.weights[i];
        store->eq.c -= store->eq.a[i];
    }
    store->eq.c = store->eq.c * 2 / TWO_TIMES_R;

    return &store->eq;
}
