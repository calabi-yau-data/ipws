#include "Global.h"
#include "LG.h"
#include "Rat.h"

FILE *inFILE, *outFILE;

int IsDigit(char c)
{
    return (('0' <= c) && (c <= '9'));
}

void RgcWeights(int narg, char *fn[]);

int main(int narg, char *fn[])
{
    inFILE = stdin;
    outFILE = stdout;

    if (narg == 1) {
        printf("\nFor help type `%s -h'\n\n", fn[0]);
        exit(0);
    }
    if (fn[1][1] == 'd')
        RgcWeights(narg, fn);
    return 0;
}

#define WDIM 800000

typedef struct {
    int d, r2, allow11; // Classification parameters
    Long x[POLY_Dmax + 1][POLY_Dmax]; // List of points that have to be allowed
                                      // by the weight system
    EqList q[POLY_Dmax]; // TODO: Precursor to weight systems. Written in
                         // ComputeQ0 and ComputeQ
    INCI qI[POLY_Dmax][EQUA_Nmax];
    int f0[POLY_Dmax]; // TODO: This is something to check for redundant points
    Equation wli[WDIM]; // Unique weight system candidates
    Long wnum; // Number of weight system candidates
    Long candnum; // Number of weight system candidates, including duplicates
    Long winum; // Number of IP weight systems
} RgcClassData;

// Compares weight systems
int RgcWeicomp(Equation w1, Equation w2, int d)
{
    /* w2-w1, i.e. pos for w1<w2,neg for w1>w2  */
    int i = d - 1;
    if (w1.c - w2.c)
        return w1.c - w2.c;
    while ((i) && (w1.a[i] == w2.a[i]))
        i--;
    return w2.a[i] - w1.a[i];
}

void RgcInsertat(Equation ww, int position, RgcClassData *X)
{
    int i;
    for (i = X->wnum - 1; i >= position; i--)
        X->wli[i + 1] = X->wli[i];
    X->wli[position] = ww;
    X->wnum++;
}

// Adds weight system wn to the sorted list X->wli, if it is basic and not
// already there.
void RgcAddweight(Equation wn, RgcClassData *X)
{
    int i, j, p, n0, n1, k;

    // Check buffer size
    if (X->wnum >= WDIM) {
        if (X->wnum > WDIM)
            return;
        else {
            X->wnum++;
            printf("WDIM too small!\n");
            fflush(0);
            return;
        }
    }

    // Skip weight systems containing a weight of 1/2
    for (i = 0; i < X->d; i++)
        if (2 * wn.a[i] == -wn.c)
            return;

    // Skip weight systems containing a weight of 1
    for (i = 0; i < X->d; i++)
        if (wn.a[i] == -wn.c)
            return;

    // Skip weight systems containing two weights with a sum of 1
    if (!X->allow11)
        for (i = 0; i < X->d - 1; i++)
            for (j = i + 1; j < X->d; j++)
                if (wn.a[i] + wn.a[j] == -wn.c)
                    return;

    // Sort weights
    X->candnum++;
    for (i = 0; i < X->d - 1; i++)
        for (p = i + 1; p < X->d; p++)
            if (wn.a[i] > wn.a[p]) {
                k = wn.a[i];
                wn.a[i] = wn.a[p];
                wn.a[p] = k;
            } /* make n0<=n1<=...<=n# */

    // Add weight system to the sorted list if it is not already there
    if (X->wnum) {
        i = RgcWeicomp(wn, X->wli[n0 = 0], X->d);
        if (!i)
            return;
        if (i > 0) {
            RgcInsertat(wn, 0, X);
            return;
        }
        i = RgcWeicomp(wn, X->wli[n1 = X->wnum - 1], X->d);
        if (!i)
            return;
        if (i < 0) {
            RgcInsertat(wn, X->wnum, X);
            return;
        }
        while (n1 > n0 + 1) {
            p = (n0 + n1) / 2;
            i = RgcWeicomp(wn, X->wli[p], X->d);
            if (!i)
                return;
            if (i > 0)
                n1 = p;
            else
                n0 = p;
        }
        RgcInsertat(wn, n1, X);
    } else
        RgcInsertat(wn, 0, X);
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
        for (i = 0; i < X->d; i++)
            printf("%d ", (int)X->q[n].e[j].a[i]);
        printf("\n");
    }
}

void PrintEquation(Equation *q, int d /*, char *c, int j*/)
{
    int i;
    printf("%d  ", (int)-q->c);
    for (i = 0; i < d; i++)
        printf("%d ", (int)q->a[i]);
    /*printf("  %s  np=%d\n", c, j);*/
}

// Tests if the last point added to X->x should really be considered.
int LastPointForbidden(int n, RgcClassData *X)
{
    int l;
    Long *y = X->x[n];
    Long ysum = 0, ymax = 0;

    assert(n < X->d);

    for (l = 0; l < X->d; l++) {
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
    if (X->r2 == 2)
        if (ymax < 2)
            return 1;

    // TODO: Why can we exclude this?
    if (X->r2 == 1)
        if (ymax < 3)
            return 1;

    // TODO: This drastically removes redundant points. How does it work?
    for (l = X->f0[n - 1]; l < X->d - 1; l++)
        if (y[l] < y[l + 1])
            return 1;

    return 0;
}

// Initializes X->q such that it represents the following weight systems:
// (r, 0, ... 0), (0, r, ... 0), ... (0, 0, ... r).
// TODO: Also does something to X->f0 and X->qI
void ComputeQ0(RgcClassData *X)
{
    int i, j;
    X->q[0].ne = X->d;
    for (i = 0; i < X->q[0].ne; i++) {
        for (j = 0; j < X->d; j++)
            X->q[0].e[i].a[j] = 0;
        if (X->r2 % 2) {
            X->q[0].e[i].a[i] = X->r2;
            X->q[0].e[i].c = -2;
        } else {
            X->q[0].e[i].a[i] = X->r2 / 2;
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
    assert(n < X->d);
    for (i = X->d - 1; (i >= X->f0[n - 1]) && (y[i] == 0); i--)
        ;
    X->f0[n] = ++i;
    qNew->ne = 0;
    for (i = 0; i < qOld->ne; i++)
        if (!(yqOld[i] = Eval_Eq_on_V(&(qOld->e[i]), y, X->d))) {
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
                        qNew->e[qNew->ne] =
                            EEV_To_Equation(&qOld->e[i], &qOld->e[j], y, X->d);
                        if (qNew->e[qNew->ne].c > 0) {
                            for (k = 0; k < X->d; k++)
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
void Cancel(Equation *q, int d)
{
    Long gcd = -q->c;
    int j;
    for (j = 0; j < d; j++)
        gcd = Fgcd(gcd, q->a[j]);
    if (gcd > 1) {
        for (j = 0; j < d; j++)
            q->a[j] /= gcd;
        q->c /= gcd;
    }
}

int ComputeAndAddAverageWeight(Equation *q, int n, RgcClassData *X)
{
    int i, j;
    EqList *el = &X->q[n];

    if (el->ne < X->d - n)
        return 0;

    q->c = -1;
    for (i = 0; i < el->ne; i++) {
        if (el->e[i].c >= 0) {
            PrintQ(n, X);
            exit(0);
        }
        q->c = -Flcm(-q->c, -el->e[i].c);
    }
    for (j = 0; j < X->d; j++) {
        q->a[j] = 0;
        for (i = 0; i < el->ne; i++)
            q->a[j] += el->e[i].a[j] * (q->c / el->e[i].c);
        if (q->a[j] <= 0) {
            assert(q->a[j] == 0);
            return 0;
        }
    }
    q->c *= el->ne;
    Cancel(q, X->d);
    RgcAddweight(*q, X);
    return 1;
}

void ComputeAndAddLastQ(RgcClassData *X)
{
    /* q[d-1] from q[d-2], x[d-1] */
    int i, d = X->d;
    Equation q;
    Equation *q0 = &X->q[d - 2].e[0], *q1 = &X->q[d - 2].e[1];
    Long *y = X->x[d - 1];
    Long yq0 = Eval_Eq_on_V(q0, y, d), yq1;
    if (LastPointForbidden(d - 1, X))
        return;
    if (!yq0)
        return;
    yq1 = Eval_Eq_on_V(q1, y, d);
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
    for (i = 0; i < d; i++)
        q.a[i] = yq0 * q1->a[i] + yq1 * q0->a[i];
    Cancel(&q, d);
    RgcAddweight(q, X);
}

void RecConstructRgcWeights(int n, RgcClassData *X)
{
    /* we have q[n-1], x[n] */
    int k, l;
    Equation q;
    Long yq[POLY_Dmax];
    Long *y = X->x[n + 1];
    if (n == 0)
        ComputeQ0(X);
    else if (LastPointForbidden(n, X))
        return;
    else
        ComputeQ(n, X);
    if (!ComputeAndAddAverageWeight(&q, n, X))
        return;
    if (n >= X->d - 1)
        return;
    /* Examine all integer points of simplex:                               */
    for (k = 0; k < X->d - 1; k++) {
        y[k] = 0;
        yq[k] = 0;
    }          /* sets k=d-1; important!    */
    y[k] = -1; /* starting point just outside                       */
    yq[k] = -q.a[k];
    while (k >= 0) {
        y[k]++;
        yq[k] += q.a[k];
        for (l = k + 1; l < X->d; l++)
            yq[l] = yq[k];
        if (n == X->d - 2) {
            assert(X->q[X->d - 2].ne == 2);
            ComputeAndAddLastQ(X);
        } else
            RecConstructRgcWeights(n + 1, X);
        for (k = X->d - 1; (k >= 0 ? (yq[k] + q.a[k] >= -q.c) : 0); k--)
            y[k] = 0;
    }
    /* sets k to the highest value where y[k] didn't exceed max value;
       resets the following max values to min values                 */
}

void AddPointToPoly(Long *y, PolyPointList *P)
{
    int i;
    for (i = 0; i < P->n; i++)
        P->x[P->np][i] = y[i];
    P->np++;
}

int WsIpCheck(Equation *q, int d)
{
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
    y[k] = -1; /* starting point just outside                       */
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
    /* sets k to the highest value where y[k] didn't exceed max value;
       resets the following max values to min values                 */
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
    k = P->np - 1;
    free(P);
    return k;
}

void RgcWeights(int narg, char *fn[])
{
    int i, j, d, n = 1, r2 = 2;
    char *c = &fn[1][2];
    RgcClassData *X = (RgcClassData *)malloc(sizeof(RgcClassData));
    if (narg > 2)
        if (c[0] == 0)
            c = fn[++n];
    if (!IsDigit(c[0])) {
        puts("-d must be followed by a number");
        exit(0);
    }
    if (POLY_Dmax < (d = atoi(c))) {
        printf("Increase POLY_Dmax to %d\n", d);
        exit(0);
    }
    if (narg > ++n) {
        if ((fn[n][0] != '-') || (fn[n][1] != 'r')) {
            printf("the second option has to be of the type -r\n");
            exit(0);
        }
        c = &fn[n][2];
        r2 = atoi(c);
    }
    X->d = d;
    X->r2 = r2;
    X->wnum = 0;
    X->winum = 0;
    X->candnum = 0;
    X->allow11 = 0;
    RecConstructRgcWeights(0, X);
    if (X->wnum <= WDIM) {
        for (i = 0; i < X->wnum; i++) {
            j = WsIpCheck(&X->wli[i], d);
            if (j) {
                PrintEquation(&X->wli[i], X->d);
                printf("  np=%d\n", j);
                X->winum++;
           /*else PrintEquation(&X->wli[i], X->d, "n");*/}
        }
    }
    printf("#ip=%ld, #cand=%ld(%ld)\n", X->winum, X->wnum, X->candnum);
}
