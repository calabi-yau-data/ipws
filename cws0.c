#include <time.h>

#include "Global.h"
#include "Rat.h"
#include "weight_system_store.h"

FILE *inFILE, *outFILE;

#define SKEW 0

void RgcWeights(void);

int main()
{
    inFILE = stdin;
    outFILE = stdout;

    RgcWeights();
    return 0;
}

typedef struct {
    int allow11;                      // Classification parameters
    Long x[POLY_Dmax + 1][POLY_Dmax]; // List of points that have to be allowed
                                      // by the weight system
    EqList q[POLY_Dmax]; // TODO: Precursor to weight systems. Written in
                         // ComputeQ0 and ComputeQ
    INCI qI[POLY_Dmax][EQUA_Nmax];
    int f0[POLY_Dmax]; // TODO: This is something to check for redundant points
    weight_system_store_t *wli;
    Long wnum;    // Number of weight system candidates
    Long candnum; // Number of weight system candidates, including duplicates
    Long winum;   // Number of IP weight systems
    size_t recursion_level_counts[POLY_Dmax];
    size_t weight_counts[POLY_Dmax];
    time_t start_time;
} RgcClassData;

void print_stats(const RgcClassData *X)
{
    printf("%7.2f: ", (time(NULL) - X->start_time) / 60.0);
    for (size_t i = 1; i < DIMENSION - 1; ++i)
        printf("%ld ", X->recursion_level_counts[i]);
    printf("-- ");

    for (size_t i = 1; i < DIMENSION; ++i)
        printf("%ld ", X->weight_counts[i]);
    printf("-- ");

    printf("%ld -- %d ", X->candnum, weight_system_store_size(X->wli));
    printf("\n");
    fflush(stdout);
}

// __attribute__ ((noinline))
static Long Eval_Eq(Equation *E, Long *V)
{
    Long p = E->c;
    for (int i = 0; i < DIMENSION; ++i)
        p += V[i] * E->a[i];
    return p;
}

// Adds weight system wn to the sorted list X->wli, if it is basic and not
// already there.
void RgcAddweight(Equation wn, RgcClassData *X)
{
    int i, j, p, k;

    // Skip weight systems containing a weight of 1/2
    for (i = 0; i < DIMENSION; i++)
        if (2 * wn.a[i] == -wn.c)
            return;

    // Skip weight systems containing a weight of 1
    for (i = 0; i < DIMENSION; i++)
        if (wn.a[i] == -wn.c)
            return;

    // Skip weight systems containing two weights with a sum of 1
    if (!X->allow11)
        for (i = 0; i < DIMENSION - 1; i++)
            for (j = i + 1; j < DIMENSION; j++)
                if (wn.a[i] + wn.a[j] == -wn.c)
                    return;

    // Sort weights
    X->candnum++;
    for (i = 0; i < DIMENSION - 1; i++)
        for (p = i + 1; p < DIMENSION; p++)
            if (wn.a[i] > wn.a[p]) {
                k = wn.a[i];
                wn.a[i] = wn.a[p];
                wn.a[p] = k;
            } /* make n0<=n1<=...<=n# */

    // Add weight system to the sorted list if it is not already there
    weight_system_store_insert(X->wli, &wn);
}

void PrintQ(int n, RgcClassData *X)
{
    int i, j;
    assert(n < POLY_Dmax);
    for (i = 0; i < n; i++)
        printf(" ");
    printf("q: ne=%d\n", X->q[n].ne);
    for (j = 0; j < X->q[n].ne; j++) {
        printf("%d  ", (int)X->q[n].e[j].c);
        for (i = 0; i < DIMENSION; i++)
            printf("%d ", (int)X->q[n].e[j].a[i]);
        printf("\n");
    }
}

void PrintEquation(const Equation *q /*, char *c, int j*/)
{
    int i;
    printf("%d  ", (int)-q->c);
    for (i = 0; i < DIMENSION; i++)
        printf("%d ", (int)q->a[i]);
    /*printf("  %s  np=%d\n", c, j);*/
}

// Tests if the last point added to X->x should really be considered.
int LastPointForbidden(int n, RgcClassData *X)
{
    int l;
    Long *y = X->x[n];
    Long ysum = 0, ymax = 0;

    // assert(n < DIMENSION);

    for (l = 0; l < DIMENSION; l++) {
        ysum += y[l];
        if (y[l] > ymax)
            ymax = y[l];
    }

    // Point leads to weight systems containing a weight of 1
    if (ysum < 2)
        return 1;

    // Point leads to weight systems containing a weight of 1/2 or two weights
    // with a sum of 1
    if (ysum == 2)
        if ((!X->allow11) || (ymax == 2))
            return 1;

    // Point does not allow positive weight systems (except if all coordinates
    // are 1)
    if (TWO_TIMES_R == 2)
        if (ymax < 2)
            return 1;

    // TODO: Why can we exclude this?
    if (TWO_TIMES_R == 1)
        if (ymax < 3)
            return 1;

    // TODO: This drastically removes redundant points. How does it work?
#if SKEW == 0
    for (l = X->f0[n - 1]; l < DIMENSION - 1; l++)
        if (y[l] < y[l + 1])
            return 1;
#endif

    return 0;
}

// Initializes X->q such that it represents the following weight systems:
// (r, 0, ... 0), (0, r, ... 0), ... (0, 0, ... r).
// TODO: Also does something to X->f0 and X->qI
void ComputeQ0(RgcClassData *X)
{
    int i, j;
    X->q[0].ne = DIMENSION;
    for (i = 0; i < X->q[0].ne; i++) {
        for (j = 0; j < DIMENSION; j++)
            X->q[0].e[i].a[j] = 0;
        if (TWO_TIMES_R % 2) {
            X->q[0].e[i].a[i] = TWO_TIMES_R;
            X->q[0].e[i].c = -2;
        } else {
            X->q[0].e[i].a[i] = TWO_TIMES_R / 2;
            X->q[0].e[i].c = -1;
        }
    }
    X->f0[0] = 0;
    X->qI[0][0] = INCI_1();
    for (i = 1; i < X->q[0].ne; i++)
        X->qI[0][i] = INCI_PN(X->qI[0][i - 1], 1);
}

int IsRedundant(INCI newINCI, INCI *qINew, int ne)
{
    int i;
    for (i = 0; i < ne; i++)
        if (INCI_LE(qINew[i], newINCI))
            return 1;
    return 0;
}

void ComputeQ(int n, RgcClassData *X)
{
    /* q[n] from q[n-1], x[n] */
    int i, j, k;
    Long *y = X->x[n];
    Long yqOld[EQUA_Nmax];
    EqList *qOld = &X->q[n - 1], *qNew = &X->q[n];
    INCI *qIOld = X->qI[n - 1], *qINew = X->qI[n];
    INCI newINCI;
    assert(n < DIMENSION);
    for (i = DIMENSION - 1; (i >= X->f0[n - 1]) && (y[i] == 0); i--)
        ;
    X->f0[n] = ++i;
    qNew->ne = 0;
    for (i = 0; i < qOld->ne; i++)
        if (!(yqOld[i] = Eval_Eq(&(qOld->e[i]), y))) {
            qNew->e[qNew->ne] = qOld->e[i];
            qINew[qNew->ne] = qIOld[i];
            (qNew->ne)++;
        }
    for (i = 0; i < qOld->ne - 1; i++)
        for (j = i + 1; j < qOld->ne; j++)
            if (yqOld[i] * yqOld[j] < 0)
                if (INCI_abs(newINCI = INCI_OR(qIOld[i], qIOld[j])) <= n + 1)
                    if (!IsRedundant(newINCI, qINew, qNew->ne)) {
                        for (k = qNew->ne - 1; k >= 0; k--)
                            if (INCI_LE(newINCI, qINew[k])) {
                                qINew[k] = qINew[qNew->ne - 1];
                                qNew->e[k] = qNew->e[qNew->ne - 1];
                                qNew->ne--;
                            }
                        assert(qNew->ne < EQUA_Nmax - 1);
                        qINew[qNew->ne] = newINCI;
                        qNew->e[qNew->ne] = EEV_To_Equation(
                            &qOld->e[i], &qOld->e[j], y, DIMENSION);
                        if (qNew->e[qNew->ne].c > 0) {
                            for (k = 0; k < DIMENSION; k++)
                                qNew->e[qNew->ne].a[k] *= -1;
                            qNew->e[qNew->ne].c *= -1;
                        }
                        qNew->ne++;
                    }
}

Long Flcm(Long a, Long b)
{
    return (a * b) / Fgcd(a, b);
}

// Normalizes q such that it does not have any common divisors
void Cancel(Equation *q)
{
    Long gcd = -q->c;
    int j;
    for (j = 0; j < DIMENSION; j++)
        gcd = Fgcd(gcd, q->a[j]);
    if (gcd > 1) {
        for (j = 0; j < DIMENSION; j++)
            q->a[j] /= gcd;
        q->c /= gcd;
    }
}

int ComputeAndAddAverageWeight(Equation *q, int n, RgcClassData *X)
{
    int i, j;
    EqList *el = &X->q[n];

    ++X->weight_counts[n];

    if (el->ne < DIMENSION - n)
        return 0;

    q->c = -1;
    for (i = 0; i < el->ne; i++) {
        if (el->e[i].c >= 0) {
            PrintQ(n, X);
            exit(0);
        }
        q->c = -Flcm(-q->c, -el->e[i].c);
    }
    for (j = 0; j < DIMENSION; j++) {
        q->a[j] = 0;
        for (i = 0; i < el->ne; i++)
            q->a[j] += el->e[i].a[j] * (q->c / el->e[i].c);

        q->a[j] += el->e[0].a[j] * (q->c / el->e[0].c) * SKEW;

        if (q->a[j] <= 0) {
            assert(q->a[j] == 0);
            return 0;
        }
    }

    q->c *= el->ne + SKEW;
    Cancel(q);
    RgcAddweight(*q, X);
    return 1;
}

void ComputeAndAddLastQ(RgcClassData *X)
{
    /* q[DIMENSION-1] from q[DIMENSION-2], x[DIMENSION-1] */
    int i;
    Equation q;
    Equation *q0 = &X->q[DIMENSION - 2].e[0], *q1 = &X->q[DIMENSION - 2].e[1];
    Long *y = X->x[DIMENSION - 1];
    Long yq0 = Eval_Eq(q0, y), yq1;

    ++X->weight_counts[DIMENSION - 1];
    if (++X->weight_counts[DIMENSION - 1] % 100000 == 0)
        print_stats(X);

    if (LastPointForbidden(DIMENSION - 1, X))
        return;
    if (!yq0)
        return;
    yq1 = Eval_Eq(q1, y);
    if (yq0 < 0) {
        if (yq1 <= 0)
            return;
        yq0 *= -1;
    } else {
        if (yq1 >= 0)
            return;
        yq1 *= -1;
    }
    q.c = yq0 * q1->c + yq1 * q0->c;
    for (i = 0; i < DIMENSION; i++)
        q.a[i] = yq0 * q1->a[i] + yq1 * q0->a[i];
    Cancel(&q);
    RgcAddweight(q, X);
}

void RecConstructRgcWeights(int n, RgcClassData *X)
{
    /* we have q[n-1], x[n] */
    int k, l;
    Equation q;
    Long yq[POLY_Dmax];
    Long *y = X->x[n + 1];

    ++X->recursion_level_counts[n];

    if (n == 0)
        ComputeQ0(X);
    else if (LastPointForbidden(n, X))
        return;
    else
        ComputeQ(n, X);

    print_stats(X);

    if (!ComputeAndAddAverageWeight(&q, n, X))
        return;
    if (n >= DIMENSION - 1)
        return;

    switch (n) {
    case 0:
        break;
    case 1:
        fprintf(stderr, "\n* ");
        break;
    default:
        fprintf(stderr, "%d:%d ", n, weight_system_store_size(X->wli));
        fflush(stderr);
        break;
    }

    /* Examine all integer points of simplex:                               */
    for (k = 0; k < DIMENSION - 1; k++) {
        y[k] = 0;
        yq[k] = 0;
    }          /* sets k=DIMENSION-1; important!    */
    y[k] = -1; /* starting point just outside                       */
    yq[k] = -q.a[k];
    while (k >= 0) {
        y[k]++;
        yq[k] += q.a[k];
        for (l = k + 1; l < DIMENSION; l++)
            yq[l] = yq[k];
        if (n == DIMENSION - 2) {
            assert(X->q[DIMENSION - 2].ne == 2);
            ComputeAndAddLastQ(X);
        } else
            RecConstructRgcWeights(n + 1, X);
        for (k = DIMENSION - 1; (k >= 0 ? (yq[k] + q.a[k] >= -q.c) : 0); k--)
            y[k] = 0;
    }
    /* sets k to the highest value where y[k] didn't exceed max value;
       resets the following max values to min values                 */
}

void AddPointToPoly(Long *y, PolyPointList *P)
{
    if (P->np >= POINT_Nmax) {
        printf("ohh no\n");
        exit(1);
    }

    int i;
    for (i = 0; i < P->n; i++)
        P->x[P->np][i] = y[i];
    P->np++;
}

int WsIpCheck(const Equation *q)
{
    int k, l;
    PolyPointList *P = (PolyPointList *)malloc(sizeof(PolyPointList));
    assert(P != NULL);
    VertexNumList V;
    EqList E;
    Long y[POLY_Dmax];
    Long yq[POLY_Dmax];
    P->n = DIMENSION;
    P->np = 0;
    for (k = 0; k < DIMENSION; k++) {
        y[k] = 0;
        yq[k] = 0;
    }
    k = DIMENSION - 1;
    AddPointToPoly(y, P);
    y[k] = -1; /* starting point just outside                       */
    yq[k] = -q->a[k];
    while (k >= 0) {
        y[k]++;
        yq[k] += q->a[k];
        for (l = k + 1; l < DIMENSION; l++)
            yq[l] = yq[k];
        if (yq[k] == -q->c)
            AddPointToPoly(y, P);
        for (k = DIMENSION - 1; (k >= 0 ? (yq[k] + q->a[k] > -q->c) : 0); k--)
            y[k] = 0;
    }
    /* sets k to the highest value where y[k] didn't exceed max value;
       resets the following max values to min values                 */
    if (P->np <= DIMENSION) {
        free(P);
        return 0;
    }
    Find_Equations(P, &V, &E);
    if (E.ne < DIMENSION) {
        free(P);
        return 0;
    }
    for (k = 0; k < DIMENSION; k++)
        y[k] = 1;
    for (k = 0; k < E.ne; k++)
        if (Eval_Eq(&(E.e[k]), y) <= 0)
            if (!E.e[k].c) {
                free(P);
                return 0;
            }
    k = P->np - 1;
    free(P);
    return k;
}

void RgcWeights(void)
{
    RgcClassData *X = (RgcClassData *)malloc(sizeof(RgcClassData));

    X->wnum = 0;
    X->winum = 0;
    X->candnum = 0;
    X->allow11 = 0;
    X->wli = weight_system_store_new();
    for (size_t i = 0; i < DIMENSION; ++i) {
        X->weight_counts[i] = 0;
        X->recursion_level_counts[i] = 0;
    }
    X->start_time = time(NULL);

    RecConstructRgcWeights(0, X);
    fprintf(stderr, "\n");

    X->wnum = weight_system_store_size(X->wli);

    const Equation *e;
    weight_system_store_begin_iteration(X->wli);
    while (e = weight_system_store_next(X->wli)) {
        int i = WsIpCheck(e);
        if (i) {
            PrintEquation(e);
            printf("  np=%d\n", i);
            X->winum++;
            /*else PrintEquation(&X->wli[i], DIMENSION, "n");*/
        }
    }
    printf("#ip=%ld, #cand=%ld(%ld)\n", X->winum, X->wnum, X->candnum);
}
