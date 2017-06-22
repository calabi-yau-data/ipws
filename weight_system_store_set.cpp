#include "weight_system_store.h"

#include <set>

class Comparator {
public:
    bool operator()(const Equation &w1, const Equation &w2) const
    {
        int i = DIMENSION - 1;
        if (w1.c - w2.c)
            return w1.c > w2.c;
        while ((i) && (w1.a[i] == w2.a[i]))
            i--;
        return w2.a[i] > w1.a[i];
    }
};

using Set = std::set<Equation, Comparator>;

struct weight_system_store {
    Set s;
};

static Set::iterator it;

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
    store->s.insert(*e);
}

int weight_system_store_size(weight_system_store_t *store)
{
    return store->s.size();
}

void weight_system_store_begin_iteration(weight_system_store_t *store)
{
    it = store->s.begin();
}

const Equation *weight_system_store_next(weight_system_store_t *store)
{
    if (it == store->s.end())
        return nullptr;
    return &*it++;
}
