#ifndef BLUBBEL
#define BLUBBEL

#include "Global.h"

#ifdef __cplusplus
extern "C" {
#endif

struct weight_system_store;
typedef struct weight_system_store weight_system_store_t;

weight_system_store_t *weight_system_store_new();
void weight_system_store_free(weight_system_store_t *store);

void weight_system_store_insert(weight_system_store_t *store,
                                const Equation *e);
int weight_system_store_size(weight_system_store_t *store);
void weight_system_store_begin_iteration(weight_system_store_t *store);
const Equation *weight_system_store_next(weight_system_store_t *store);

#ifdef __cplusplus
}
#endif

#endif
