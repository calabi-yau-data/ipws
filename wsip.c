#include "palp/Global.h"

FILE *inFILE, *outFILE;

static PolyPointList *P;

static int AddPointToPoly(Long *y, PolyPointList *P)
{
    if (P->np >= POINT_Nmax)
        return 0;

    int i;
    for (i = 0; i < P->n; i++)
        P->x[P->np][i] = y[i];
    P->np++;

    return 1;
}

static int WsIpCheck(Equation *q, int d)
{
    int k, l;

    if (!P) {
        P = (PolyPointList *)malloc(sizeof(PolyPointList));
        if (!P)
            return -1;
    }

    VertexNumList V;
    EqList E;
    Long y[POLY_Dmax];
    Long yq[POLY_Dmax];
    P->n = d;
    P->np = 0;
    for (k = 0; k < d; k++) {
        y[k] = 0;
        yq[k] = 0;
    }
    k = d - 1;
    if (!AddPointToPoly(y, P))
        return -1;
    y[k] = -1; // starting point just outside
    yq[k] = -q->a[k];
    while (k >= 0) {
        y[k]++;
        yq[k] += q->a[k];
        for (l = k + 1; l < d; l++)
            yq[l] = yq[k];
        if (yq[k] == -q->c)
            if (!AddPointToPoly(y, P))
                return -1;

        for (k = d - 1; (k >= 0 ? (yq[k] + q->a[k] > -q->c) : 0); k--)
            y[k] = 0;
    }
    if (P->np <= d)
        return 0;
    Find_Equations(P, &V, &E);
    if (E.ne < d)
        return 0;
    for (k = 0; k < d; k++)
        y[k] = 1;
    for (k = 0; k < E.ne; k++)
        if (Eval_Eq_on_V(&(E.e[k]), y, d) <= 0)
            if (!E.e[k].c)
                return 0;

    return P->np != 1 ? 1 : 0;
}

__attribute__((visibility("default")))
int weight_system_has_ip(const int *weights, int divisor, int dim)
{
    if (dim > POLY_Dmax)
        return -1;

    Equation eq;

    for (int i = 0; i < dim; ++i)
        eq.a[i] = weights[i];
    eq.c = -divisor;

    return WsIpCheck(&eq, dim);
}
