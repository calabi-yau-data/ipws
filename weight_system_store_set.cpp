#include "weight_system_store.h"

#include <set>

class WeightSystem {
    int32_t weights[dim];

public:
    bool operator<(const WeightSystem &rhs) const
    {
        for (size_t i = dim; i--;)
            if (weights[i] != rhs.weights[i])
                return weights[i] < rhs.weights[i];
        return false;
    }

    WeightSystem(const Hyperplane &eq)
    {
        for (size_t i = 0; i < dim; ++i) {
            assert(eq.a[i] >= 0 && eq.a[i] <= INT32_MAX);
            weights[i] = static_cast<int32_t>(eq.a[i]);
        }
    }

    operator Hyperplane() const
    {
        Hyperplane eq;

        eq.c = 0;
        for (size_t i = 0; i < dim; ++i) {
            eq.a[i] = weights[i];
            eq.c -= eq.a[i];
        }
        eq.c = eq.c * 2 / two_times_r;

        return eq;
    }
};

using Set = std::set<WeightSystem>;

struct weight_system_store {
    Set s;
    Set::iterator it;
    Hyperplane eq;
};

weight_system_store_t *weight_system_store_new()
{
    return new weight_system_store_t{};
}

void weight_system_store_free(weight_system_store_t *store)
{
    delete store;
}

void weight_system_store_insert(weight_system_store_t *store,
                                const Hyperplane *e)
{
    store->s.insert(WeightSystem{*e});
}

int weight_system_store_size(weight_system_store_t *store)
{
    return static_cast<int>(store->s.size());
}

void weight_system_store_begin_iteration(weight_system_store_t *store)
{
    store->it = store->s.begin();
}

const Hyperplane *weight_system_store_next(weight_system_store_t *store)
{
    if (store->it == store->s.end())
        return nullptr;

    const WeightSystem &ws = *store->it++;
    store->eq = ws;

    return &store->eq;
}
