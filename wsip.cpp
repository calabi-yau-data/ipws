#include "global.h"

int Find_Equations(PolyPointList *P, VertexNumList *VNL, EqList *EL);

Long Eval_Eq_on_V(Equation *E, Long *V, int i)
{
    Long p = E->c;
    while (i--)
        p += V[i] * E->a[i];
    return p;
}

namespace {

#define MAX_BAD_EQ (dim > 5) /* previously 6; needed for nef !? */
#define CEQ_Nmax EQUA_Nmax

Long Egcd(Long A0, Long A1, Long *Vout0, Long *Vout1)
{
    Long V0 = A0, V1 = A1, A2, X0 = 1, X1 = 0, X2 = 0;
    while ((A2 = A0 % A1)) {
        X2 = X0 - X1 * (A0 / A1);
        A0 = A1;
        A1 = A2;
        X0 = X1;
        X1 = X2;
    }
    *Vout0 = X1, *Vout1 = (A1 - (V0)*X1) / (V1);
    return A1;
}
Long NNgcd(Long a, Long b) /* NonNegative gcd handling 0 */
{
    a = (a < 0) ? -a : a;
    b = (b < 0) ? -b : b;
    if (!b)
        return a;
    while (a %= b)
        if (!(b %= a))
            return a;
    return b;
}
LLong LNNgcd(LLong a, LLong b) /* NonNegative gcd handling 0 */
{
    a = (a < 0) ? -a : a;
    b = (b < 0) ? -b : b;
    if (!b)
        return a;
    while (a %= b)
        if (!(b %= a))
            return a;
    return b;
}

typedef struct {
    int ne;
    Equation e[CEQ_Nmax];
} CEqList;

Long RoundQ(Long N, Long D)
{
    Long F;
    if (D < 0) {
        D = -D;
        N = -N;
    }
    F = N / D;
    return F + (2 * (N - F * D)) / D;
}

Long W_to_GLZ(Long *W, int *d, Long **GLZ)
{
    int i, j;
    Long G, *E = *GLZ, *B = GLZ[1];
    for (i = 0; i < *d; i++)
        assert(W[i] != 0);
    for (i = 1; i < *d; i++)
        for (j = 0; j < *d; j++)
            GLZ[i][j] = 0;
    G = Egcd(W[0], W[1], &E[0], &E[1]);
    B[0] = -W[1] / G;
    B[1] = W[0] / G;
    for (i = 2; i < *d; i++) {
        Long a, b, g = Egcd(G, W[i], &a, &b);
        B = GLZ[i];
        B[i] = G / g;
        G = W[i] / g;
        for (j = 0; j < i; j++)
            B[j] = -E[j] * G; /* B=Base-Line */
        for (j = 0; j < i; j++)
            E[j] *= a;
        E[j] = b;                   /* next Egcd */
        for (j = i - 1; 0 < j; j--) /* I M P R O V E M E N T */
        {
            int n;
            Long *Y = GLZ[j], rB = RoundQ(B[j], Y[j]), rE = RoundQ(E[j], Y[j]);
            /*  rB=CeilQ(B[j],Y[j]), rE=CeilQ(E[j],Y[j]); */
            /*  printf(" [%d/%d -> %d] ",(int)B[j],(int)Y[j],(int)rB);
                printf(" [%d/%d -> %d] ",(int)E[j],(int)Y[j],(int)rE); 	*/
            for (n = 0; n <= j; n++) {
                B[n] -= rB * Y[n];
                E[n] -= rE * Y[n];
            }
        }
        G = g;
    }
    return G;
}

int INCI_abs(INCI X)
{
    int abs = 0;
    while (!INCI_EQ_0(X)) {
        abs += INCI_M2(X);
        X = INCI_D2(X);
    }
    return abs;
}

INCI Eq_To_INCI(Equation *_Eq, PolyPointList *_P, VertexNumList *_V)
{
    int j;
    INCI X = INCI_0();
    for (j = 0; j < _V->nv; j++)
        X = INCI_PN(X, Eval_Eq_on_V(_Eq, _P->x[_V->v[j]], dim));
    return X;
}

/*  ======================================================================  */
/*  ==========	     			   	  	  	==========  */
/*  ========== G E N E R A L   P U R P O S E   R O U T I N E S  ==========  */
/*  ==========						  	==========  */
/*  ======================================================================  */

#define LLong_EEV (1) /* 1 @ [4662 4 20 333 422 1554 2329] */
#define TEST_EEV (0)  /* compare Long to LLong EEV */

Equation EEV_To_Equation(Equation *_E1, Equation *_E2, Long *_V, int n)
{
    /* Calculate the equation spanned by _V and the intersection of _E1, _E2  */
    int i;
    Long l, m, g;
    Equation Eq;
    l = Eval_Eq_on_V(_E2, _V, n);
    m = Eval_Eq_on_V(_E1, _V, n);
    g = NNgcd(l, m);
    assert(g);
    l /= g;
    m /= g;
#if ((!(LLong_EEV)) || (TEST_EEV)) /* Long version */
    for (i = 0; i < n; i++)
        Eq.a[i] = l * _E1->a[i] - m * _E2->a[i];
    {
        int gcd = Eq.c = l * _E1->c - m * _E2->c;
        for (i = 0; i < n; i++)
            gcd = NNgcd(gcd, Eq.a[i]);
        assert(gcd);
        if (gcd != 1) {
            for (i = 0; i < n; i++)
                Eq.a[i] /= gcd;
            Eq.c /= gcd;
        }
    }
#endif
#if ((LLong_EEV) || (TEST_EEV)) /* LLong version */
    {
        LLong A[dim], C, G;
        for (i = 0; i < n; i++)
            A[i] = ((LLong)l) * ((LLong)_E1->a[i]) -
                   ((LLong)m) * ((LLong)_E2->a[i]);
        G = C = ((LLong)l) * ((LLong)_E1->c) - ((LLong)m) * ((LLong)_E2->c);
        for (i = 0; i < n; i++)
            G = LNNgcd(G, A[i]);
        assert(G);
        if (G != 1) {
            C /= G;
            for (i = 0; i < n; i++)
                A[i] /= G;
        }
        Eq.c = C;
        for (i = 0; i < n; i++)
            Eq.a[i] = A[i];
    }
#endif
    return Eq;
}

int Vec_Greater_Than(Long *X, Long *Y, int i)
{ /* return 1 iff `X > Y' */
    while (i--) {
        if (X[i] > Y[i])
            return 1;
        if (X[i] < Y[i])
            return 0;
    }
    puts("Identical points in Vec_Greater_Than !!");
    exit(0);
    return 0;
}

int IsGoodCEq(Equation *_E, PolyPointList *_P, VertexNumList *_V)
{
    int i = _V->nv;
    Long s;
    while (!(s = Eval_Eq_on_V(_E, _P->x[_V->v[--i]], dim)))
        ;
    if (s < 0) {
        int j = dim;
        while (j--)
            _E->a[j] = -_E->a[j];
        _E->c = -_E->c;
    }
    while (i)
        if (Eval_Eq_on_V(_E, _P->x[_V->v[--i]], dim) < 0)
            return 0;
    return 1;
}

int Search_New_Vertex(Equation *_E, PolyPointList *_P)
{
    int i, v = 0;
    Long *X = _P->x[0], x = Eval_Eq_on_V(_E, X, dim);
    for (i = 1; i < _P->np; i++) {
        Long *Y = _P->x[i], y = Eval_Eq_on_V(_E, Y, dim);
        if (y > x)
            continue;
        if (y == x)
            if (Vec_Greater_Than(X, Y, dim))
                continue;
        v = i;
        X = Y;
        x = y;
    }
    return v;
}

/*  ======================================================================  */
/*  ==========		     			  		==========  */
/*  ==========          S T A R T -- S I M P L E X              ==========  */
/*  ==========		     			  		==========  */
/*  ======================================================================  */

/*	return 0 <=> max.dim., E.ne==P.n+1, made Simplex of Vertices of P;  *
 *      return (P.n-E.ne) == codim. > 0  <=> E.ne defining equations on E;  */

/*  #define  NEW_START_SIMPLEX  (1)     1 @ [16644 1 38 439 2315 5548 8303] */

#define VERT_WITH_MAX_DISTANCE (0) /* 0 @ [1845 2 15 97 247 610 874]    */
#define LONG_EQ_FIRST (0)          /* 0 @ [3425 2 7 137 429 1141 1709]  */
#define TEST_GLZ_EQ (0)            /* trace StartSimplex EQs  */

Long VZ_to_Base(Long *V, int *d, Long M[dim][dim]) /* 0 iff V=0 */
{
    int p[dim], i, j, J = 0;
    Long g = 0, W[dim], *G[dim];
    for (i = 0; i < *d; i++)
        if (V[i]) {
            W[J] = V[i];
            G[J] = M[i];
            p[J++] = i;
        } else
            for (j = 0; j < *d; j++)
                M[i][j] = (i == j);
    if (J)
        if (p[0]) {
            G[0] = M[0];
            for (j = 0; j < *d; j++)
                M[p[0]][j] = (j == 0);
        }
    if (J > 1)
        g = W_to_GLZ(W, &J, G);
    else if (J) {
        g = *W;
        M[0][0] = 0;
        M[0][p[0]] = 1;
    }
    if (J > 1) {
        for (i = 0; i < J; i++) {
            int I = J;
            for (j = *d - 1; j >= 0; j--)
                G[i][j] = (V[j]) ? G[i][--I] : 0;
            assert(I == 0);
        }
    }
    return g;
}

int OrthBase_red_by_V(Long *V, Long A[][dim], int *r,
                      Long B[][dim])
{
    int i, j, k;
    Long W[dim], G[dim][dim];
    for (i = 0; i < *r; i++) {
        int j;
        W[i] = 0;
        for (j = 0; j < dim; j++)
            W[i] += A[i][j] * V[j];
    }
    assert(VZ_to_Base(W, r, G));
    for (i = 0; i < *r - 1; i++)
        for (k = 0; k < dim; k++) {
            B[i][k] = 0;
            for (j = 0; j < *r; j++)
                B[i][k] += G[i + 1][j] * A[j][k];
        }
    return (*r)--;
}

int New_Start_Vertex(Long *V0, Long *Ea, PolyPointList *P, int *v) /* P.x[v] */
{
    Equation E;
    int i, n = 0, p = 0;
    Long d, dn = 0, dp = 0, *Xn = P->x[0], *Xp = Xn;
    for (i = 0; i < dim; i++)
        E.a[i] = Ea[i];
    E.c = 0;
    E.c = -Eval_Eq_on_V(&E, V0, dim);
    d = Eval_Eq_on_V(&E, P->x[0], dim);
    if (d > 0)
        dp = d;
    if (d < 0)
        dn = d;
    for (i = 1; i < P->np; i++) {
        d = Eval_Eq_on_V(&E, P->x[i], dim);
        if (d == 0)
            continue;
        if (d == dp)
            if (Vec_Greater_Than(P->x[i], Xp, dim))
                Xp = P->x[p = i];
        if (d > dp) {
            dp = d;
            Xp = P->x[p = i];
        }
        if (d == dn)
            if (Vec_Greater_Than(P->x[i], Xn, dim))
                Xn = P->x[n = i];
        if (d < dn) {
            dn = d;
            Xn = P->x[n = i];
        }
    }
    if (dp)
        if (dn) /* points on both sides */
#if (VERT_WITH_MAX_DISTANCE)
        {
            if (dp + dn > 0)
                *v = p;
            else
                *v = n;
        }
#else
        {
            if (dp + dn > 0)
                *v = n;
            else
                *v = p;
        }
#endif
        else
            *v = p; /* d >=0 */
    else if (dn)
        *v = n; /* d <=0 */
    else
        return 0;
    /*	for(i=0;i<dim;i++) printf(" %d ",Xp[i]); printf(" = Xp  Xn =");
            for(i=0;i<dim;i++) printf(" %d ",Xn[i]); printf(" \n");
    */
    return 1;
}

/*   =========	   GLZ_Start_Simplex()  =>  return codimension 	  =======   */
int GLZ_Start_Simplex(PolyPointList *_P, VertexNumList *_V, CEqList *_C)
{
    int i, x = 0, y = 0, *VN = _V->v, r = dim, b[dim];
    Long *X = _P->x[x], *Y = _P->x[y], XX = 0, YY = 0,
         B[(dim * (dim + 1)) / 2][dim], W[dim];
    if (_P->np < 2) {
        for (x = 0; x < dim; x++)
            for (y = 0; y < dim; y++)
                _C->e[x].a[y] = (x == y);
        assert(_P->np > 0);
        for (x = 0; x < dim; x++)
            _C->e[x].c = -_P->x[0][x];
        return _C->ne = dim;
    }
    for (i = 1; i < _P->np; i++) {
        Long *Z = _P->x[i];
        if (Vec_Greater_Than(X, Z, dim))
            X = _P->x[x = i]; /* (x_n)-max: VN[0] */
        if (Vec_Greater_Than(Z, Y, dim))
            Y = _P->x[y = i]; /* (x_n)-min: VN[1] */
    }
    assert(x != y); /* at this point I need two different vertices */
    for (i = 0; i < dim; i++) {
        Long Xi = (X[i] > 0) ? X[i] : -X[i], Yi = (Y[i] > 0) ? Y[i] : -Y[i];
        if (Xi > XX)
            XX = Xi;
        if (Yi > YY)
            YY = Yi;
    }
    if (YY < XX) {
        VN[0] = y;
        VN[1] = x;
    } else {
        VN[0] = x;
        VN[1] = y;
    }
    _V->nv = 2;
    y = VN[1];
    X = _P->x[VN[0]];
    Y = _P->x[VN[1]];
    for (i = 0; i < dim; i++)
        b[i] = (i * (2 * dim - i + 1)) / 2;
    for (x = 0; x < dim; x++)
        for (i = 0; i < dim; i++)
            B[x][i] = (x == i); /* b[i+1]-b[i]=d-i */
    for (x = 1; x < dim; x++) {
        for (i = 0; i < dim; i++)
            W[i] = _P->x[y][i] - X[i];
        OrthBase_red_by_V(W, &B[b[x - 1]], &r, &B[b[x]]);
        for (i = 0; i < r; i++)
#if (LONG_EQ_FIRST)
            if (New_Start_Vertex(X, B[b[x] + r - i - 1], _P, &y))
                break;
#else
            if (New_Start_Vertex(X, B[b[x] + i], _P, &y))
                break;
#endif
        if (i == r)
            break;
        _V->v[_V->nv++] = y; /* x = dim(span) < d */
    }
    if (x < dim) {
        for (y = 0; y < r; y++) {
            Equation *E = &_C->e[y];
            Long *Z = B[b[x] + y];
            E->c = 0;
            for (i = 0; i < dim; i++)
                E->a[i] = Z[i];
            E->c = -Eval_Eq_on_V(E, X, dim);
        }
        return _C->ne = r;
    } else {
        Equation *E = _C->e;
        Long *Z = B[b[dim - 1]];
        E->c = 0;
        _C->ne = 2;
        for (i = 0; i < dim; i++)
            E->a[i] = Z[i];
        E->c = -Eval_Eq_on_V(E, X, dim);
        if (Eval_Eq_on_V(E, _P->x[_V->v[dim]], dim) < 0) {
            for (i = 0; i < dim; i++)
                E->a[i] = -Z[i];
            E->c *= -1;
        }
        X = _P->x[_V->v[r = dim]];
        for (x = 1; x < dim; x++) /* now the 2nd equation */
        {
            Y = _P->x[_V->v[x - 1]];
            for (i = 0; i < dim; i++)
                W[i] = X[i] - Y[i];
            OrthBase_red_by_V(W, &B[b[x - 1]], &r, &B[b[x]]);
        }
        E = &_C->e[1];
        E->c = 0;
        for (i = 0; i < dim; i++)
            E->a[i] = Z[i];
        E->c = -Eval_Eq_on_V(E, X, dim);
        assert(XX = Eval_Eq_on_V(E, _P->x[_V->v[dim - 1]], dim));
        if (XX < 0) {
            for (i = 0; i < dim; i++)
                E->a[i] = -Z[i];
            E->c *= -1;
        }
        for (x = dim - 2; x >= 0; x--) /* omit vertex #x */
        {
            r = dim - x;
            for (y = x + 1; y < dim; y++) {
                Y = _P->x[_V->v[y]];
                for (i = 0; i < dim; i++)
                    W[i] = X[i] - Y[i];
                OrthBase_red_by_V(W, &B[b[y - 1]], &r, &B[b[y]]);
            }
            E = &_C->e[(_C->ne)++];
            E->c = 0;
            for (i = 0; i < dim; i++)
                E->a[i] = Z[i];
            E->c = -Eval_Eq_on_V(E, X, dim);
            assert(XX = Eval_Eq_on_V(E, _P->x[_V->v[x]], dim));
            if (XX < 0) {
                for (i = 0; i < dim; i++)
                    E->a[i] = -Z[i];
                E->c *= -1;
            }
        }
    }
    assert(dim + 1 == _C->ne);
    for (x = 0; x < _C->ne; x++)
        for (i = 0; i <= dim; i++)
            assert((x == i) ==
                   (0 != Eval_Eq_on_V(&_C->e[x], _P->x[_V->v[dim - i]], dim)));
    return 0;
}

/*  ======================================================================  */
/*  ==========		     			  		==========  */
/*  ==========	     P O L Y H E D R O N   A N A L Y S I S  	==========  */
/*  ==========							==========  */
/*  ======================================================================  */

void Make_New_CEqs(PolyPointList *_P, VertexNumList *_V, CEqList *_C,
                   EqList *_F, INCI *CEq_I, INCI *F_I)
{
    int i, j, Old_C_ne = _C->ne;
    static CEqList Bad_C;
    static INCI Bad_C_I[CEQ_Nmax];

    Bad_C.ne = _C->ne = 0;
    for (i = 0; i < Old_C_ne; i++) {
        Long dist = Eval_Eq_on_V(&_C->e[i], _P->x[_V->v[_V->nv - 1]], dim);
        CEq_I[i] = INCI_PN(CEq_I[i], dist);
        if (dist < 0) {
            Bad_C.e[Bad_C.ne] = _C->e[i];
            Bad_C_I[Bad_C.ne++] = CEq_I[i];
        } else {
            _C->e[_C->ne] = _C->e[i];
            CEq_I[_C->ne++] = CEq_I[i];
        }
    }

    Old_C_ne = _C->ne;
    for (i = 0; i < _F->ne; i++)
        F_I[i] = INCI_PN(
            F_I[i], Eval_Eq_on_V(&_F->e[i], _P->x[_V->v[_V->nv - 1]], dim));
    for (j = 0; j < _F->ne; j++)
        if (!INCI_M2(F_I[j]))
            for (i = 0; i < Bad_C.ne; i++) {
                INCI New_Face = INCI_AND(Bad_C_I[i], F_I[j]);
                int k;
                if (INCI_abs(New_Face) < dim - 1)
                    continue;
                for (k = 0; k < Bad_C.ne; k++)
                    if (INCI_LE(New_Face, Bad_C_I[k]))
                        if (k != i)
                            break;
                if (k != Bad_C.ne)
                    continue;
                for (k = 0; k < Old_C_ne; k++)
                    if (INCI_LE(New_Face, CEq_I[k]))
                        break;
                if (k != Old_C_ne)
                    continue;
                for (k = 0; k < _F->ne; k++)
                    if (INCI_LE(New_Face, F_I[k]))
                        if (k != j)
                            break;
                if (k != _F->ne)
                    continue;
                assert(_C->ne < CEQ_Nmax);
                CEq_I[_C->ne] = INCI_PN(INCI_D2(New_Face), 0);
                _C->e[_C->ne] =
                    EEV_To_Equation(&(Bad_C.e[i]), &(_F->e[j]),
                                    _P->x[_V->v[_V->nv - 1]], dim);
                assert(IsGoodCEq(&(_C->e[_C->ne++]), _P, _V));
            }
    for (j = 0; j < Old_C_ne; j++)
        if (!INCI_M2(CEq_I[j]))
            for (i = Bad_C.ne - 1; i >= 0; i--) {
                INCI New_Face = INCI_AND(Bad_C_I[i], CEq_I[j]);
                int k;
                if (INCI_abs(New_Face) < dim - 1)
                    continue;
                for (k = 0; k < Bad_C.ne; k++)
                    if (INCI_LE(New_Face, Bad_C_I[k]))
                        if (k != i)
                            break;
                if (k != Bad_C.ne)
                    continue;
                for (k = 0; k < Old_C_ne; k++)
                    if (INCI_LE(New_Face, CEq_I[k]))
                        if (k != j)
                            break;
                if (k != Old_C_ne)
                    continue;
                for (k = 0; k < _F->ne; k++)
                    if (INCI_LE(New_Face, F_I[k]))
                        break;
                if (k != _F->ne)
                    continue;
                assert(_C->ne < CEQ_Nmax);
                CEq_I[_C->ne] = INCI_PN(INCI_D2(New_Face), 0);
                _C->e[_C->ne] =
                    EEV_To_Equation(&(Bad_C.e[i]), &(_C->e[j]),
                                    _P->x[_V->v[_V->nv - 1]], dim);
                assert(IsGoodCEq(&(_C->e[_C->ne++]), _P, _V));
            }
}

#if MAX_BAD_EQ

int INCI_lex_GT(INCI *x, INCI *y)
{
    return (*x > *y) ? 1 : 0;
}

int FE_Search_Bad_Eq(CEqList *_C, EqList *_F, INCI *CEq_I, INCI *F_I,
                     PolyPointList *_P, int *_IP)
{ /* return 0 :: no bad eq. */
    while (_C->ne--) {
        int j, M = _C->ne; /* INCI_LmR INCI_lex_GT */
        for (j = 0; j < _C->ne; j++)
            if (INCI_lex_GT(&CEq_I[j], &CEq_I[M]))
                M = j;
        for (j = 0; j < _P->np; j++)
            if (Eval_Eq_on_V(&(_C->e[M]), _P->x[j], dim) < 0) {
                INCI AI = CEq_I[M];
                Equation AE = _C->e[M];
                CEq_I[M] = CEq_I[_C->ne];
                _C->e[M] = _C->e[_C->ne];
                CEq_I[_C->ne] = AI;
                _C->e[_C->ne] = AE;
                return ++_C->ne;
            }
        if (_C->e[M].c < 1)
            *_IP = 0;
        assert(_F->ne < EQUA_Nmax);
        /* printf("#Feq=%d  #Ceq=%d\n",_F->ne,_C->ne); fflush(stdout); */
        _F->e[_F->ne] = _C->e[M];
        F_I[_F->ne++] = CEq_I[M];
        if (M < _C->ne) {
            _C->e[M] = _C->e[_C->ne];
            CEq_I[M] = CEq_I[_C->ne];
        }
    }
    return 0;
}

#else

int FE_Search_Bad_Eq(CEqList *_C, EqList *_F, INCI *CEq_I, INCI *F_I,
                     PolyPointList *_P, int *_IP)
{ /* return 0 :: no bad eq. */
    while (_C->ne--) {
        int j;
        for (j = 0; j < _P->np; j++)
            if (Eval_Eq_on_V(&(_C->e[_C->ne]), _P->x[j], dim) < 0)
                return ++_C->ne;
        if (_C->e[_C->ne].c < 1)
            *_IP = 0;
        assert(_F->ne < EQUA_Nmax);
        _F->e[_F->ne] = _C->e[_C->ne];
        F_I[_F->ne++] = CEq_I[_C->ne];
    }
    return 0;
}

#endif

int Finish_Find_Equations(PolyPointList *_P, VertexNumList *_V, EqList *_F,
                          CEqList *_CEq, INCI *F_I, INCI *CEq_I)
{
    int IP = 1;
    while (0 <= _CEq->ne)
        if (FE_Search_Bad_Eq(_CEq, _F, CEq_I, F_I, _P, &IP)) {
            assert(_V->nv < VERT_Nmax);
            _V->v[_V->nv++] = Search_New_Vertex(&(_CEq->e[_CEq->ne - 1]), _P);
            Make_New_CEqs(_P, _V, _CEq, _F, CEq_I, F_I);
        }
    return IP;
}
}

int Find_Equations(PolyPointList *_P, VertexNumList *_V, EqList *_F)
{
    /* return: IP, finds Vertices and Equations for _P even if not IP */
    int i;
    CEqList *CEq = (CEqList *)malloc(sizeof(CEqList));
    INCI *CEq_I = (INCI *)malloc(sizeof(INCI) * CEQ_Nmax);
    INCI *F_I = (INCI *)malloc(sizeof(INCI) * EQUA_Nmax);
    CEq->ne = 0;
    if ((CEq == NULL) || (CEq_I == NULL) || (F_I == NULL)) {
        printf("Allocation failure in Find_Equations\n");
        exit(0);
    }
    if (GLZ_Start_Simplex(_P, _V, CEq)) {
        _F->ne = CEq->ne;
        for (i = 0; i < _F->ne; i++)
            _F->e[i] = CEq->e[i];
        free(CEq);
        free(CEq_I);
        free(F_I);
        return 0;
    }
    _F->ne = 0;
    for (i = 0; i < CEq->ne; i++)
        if (INCI_abs(CEq_I[i] = Eq_To_INCI(&(CEq->e[i]), _P, _V)) < dim) {
            fprintf(stderr, "Bad CEq in Find_Equations");
            exit(0);
        }
    i = Finish_Find_Equations(_P, _V, _F, CEq, F_I, CEq_I);
    free(CEq);
    free(CEq_I);
    free(F_I);
    return i;
}
