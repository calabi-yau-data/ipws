#include "Global.h"
#include "WsIp.h"

FILE *inFILE, *outFILE;

static void AddPointToPoly(Long *y, PolyPointList *P) {
    assert(P->np < POINT_Nmax);
    int i;
    for (i = 0; i < P->n; i++)
        P->x[P->np][i] = y[i];
    P->np++;
}

static bool WsIpCheck(Equation *q, int d) {
    int k, l;
    PolyPointList *P = (PolyPointList *)malloc(sizeof(PolyPointList));
    assert(P != NULL);
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
    AddPointToPoly(y, P);
    y[k] = -1; // starting point just outside
    yq[k] = -q->a[k];
    while (k >= 0) {
        y[k]++;
        yq[k] += q->a[k];
        for (l = k + 1; l < d; l++)
            yq[l] = yq[k];
        if (yq[k] == -q->c)
            AddPointToPoly(y, P);
        for (k = d - 1; (k >= 0 ? (yq[k] + q->a[k] > -q->c) : 0); k--)
            y[k] = 0;
    }
    if (P->np <= d) {
        free(P);
        return 0;
    }
    Find_Equations(P, &V, &E);
    if (E.ne < d) {
        free(P);
        return 0;
    }
    for (k = 0; k < d; k++)
        y[k] = 1;
    for (k = 0; k < E.ne; k++)
        if (Eval_Eq_on_V(&(E.e[k]), y, d) <= 0)
            if (!E.e[k].c) {
                free(P);
                return 0;
            }

    bool ret = P->np != 1;
    free(P);
    return ret;
}

JNIEXPORT jint JNICALL
Java_WsIp_hasIp(JNIEnv *env, jclass UNUSED(cl), jintArray weightsArray,
                jint divisor) {

    int d = (*env)->GetArrayLength(env, weightsArray);
    int *weights = (*env)->GetIntArrayElements(env, weightsArray, 0);

    if (d > POLY_Dmax)
        return -1;

    Equation eq;

    for (int i = 0; i < d; ++i)
        eq.a[i] = weights[i];
    eq.c = - divisor;

    return WsIpCheck(&eq, d);
}
