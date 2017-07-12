#ifndef GLOBAL_H
#define GLOBAL_H

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "config.h"

/*
These are include files that should exist in your C library.
*/

/*  ============	basic choice of PARAMETERS	      ============  */

/*
For reflexive polytopes in 4 or less dimensions, everything should work with
Long set to 32-bit-integer and LLong set to 64 bits.
Many applications will even work with LLong at 32 bits.
For higher dimensional or complicated non-reflexive polytopes it may be
necessary to set even Long to 64 bits.
*/

#define POLY_Dmax DIMENSION

/*
POLY_Dmax should be set to the dimension of the polytopes that are analysed.
While the programs still work if POLY_Dmax is set to a higher value, they may
be considerably slowed down.
*/

#if (POLY_Dmax <= 3)
#define POINT_Nmax 40 /* max number of points	    */
#define VERT_Nmax 16  /* max number of vertices   */
#define FACE_Nmax 30  /* max number of faces      */
#define SYM_Nmax 88   /* cube: 2^D*D! plus extra  */

#elif (POLY_Dmax == 4)
#define POINT_Nmax 700 /* max number of points	    */
#define VERT_Nmax 64   /* max number of vertices   */
#define FACE_Nmax 824  /* max number of faces      */
#define SYM_Nmax 1200

#else
#define POINT_Nmax 2000000
#define VERT_Nmax 64    /* !! use optimal value !!  */
#define FACE_Nmax 10000 /* max number of faces      */
#define SYM_Nmax 46080  /* symmetry (P_1)^6: 2^6*6! */
#define EQUA_Nmax 1280  /* up to 20000 without alloc */
#endif

#ifndef EQUA_Nmax /* default setting */
#define EQUA_Nmax VERT_Nmax
#endif
/*
POINT_Nmax, VERT_Nmax and FACE_Nmax denote the maximal numbers of points,
vertices and faces, respectively.
SYM_Nmax is the maximal number of symmetries of a polytope, i.e. the order of
the finite subgroup S of the group GL(n,Z) of lattice automorphisms that leaves
a polytope invariant.
EQUA_Nmax denotes the maximal number of facets (given by equations) of a
polytope. By duality this is just the number of vertices of the dual polytope,
so it makes sense to have the default setting EQUA_Nmax = VERT_Nmax.
In applications not related to reflexive polytopes or in large dimensions a
larger value may be useful. While CPU-time is almost independent of EQUA_Nmax,
it strongly depends on VERT_Nmax/32 (i.e. use 32, 64, 96, ...).
Our settings for dimensions less than or equal to 4 are such that they work
for any reflexive polytope.
*/

#define AMBI_Dmax (5 * POLY_Dmax) /* default setting */
/*
If a polytope is determined by a combined weight system it is first realised
by an embeddeding in an ambient space of dimension (Poly-dim + number of
weight systems). AMBI_Dmax is the maximal dimension of this ambient space.
*/

#define FIB_Nmax 3000 /*NOW: 27/5/11 default setting*/
/*
Given a polytope P* it is possible to analyze the IP simplices among its
points. These simplices are given in terms of weight relations among points
of P*. FIB_Nmax is the maximal number of allowed relations.
*/

#define CD2F_Nmax FACE_Nmax
/*
Max number of codimension 2 faces.
*/

#define GL_Long Long
/*
Uses W_to_GLZ like in Rat.c
*/

extern FILE *inFILE, *outFILE;
/*
Ascii-files for input and output. If not given in the parameter list they
default to stdin and stdout, respectively.
*/

/*  ==========         Global typedefs           		==========  */

typedef struct {
    int n, np;
    Long x[POINT_Nmax][POLY_Dmax];
} PolyPointList;
/*
A list (not necessarily complete) of lattice points of a polytope.
P.x[i][j] is the j'th coordinate of the i'th lattice point.
P.n is the dimension of the polytope and P.np the number of points in the list.
*/

typedef struct {
    int v[VERT_Nmax];
    int nv;
} VertexNumList;
/*
The list of vertices of a polytope, referring to some PolyPointList P.
The j'th coordinate of the i'th vertex is then given by P.x[V.v[i]][j].
V.nv is the number of vertices of P.
*/

typedef struct {
    Long a[POLY_Dmax], c;
} Equation;
/*
This structure determines an equation of the type ax+c=0, explicitly:
sum_{i=1}^n E.a[i] x_i + E.c = 0.
*/

typedef struct {
    int ne;
    Equation e[EQUA_Nmax];
} EqList;
/*
A list of equations; EL.ne is the number of equations in the list.
*/

typedef struct {
    EqList B;
    Long W[AMBI_Dmax][AMBI_Dmax], d[AMBI_Dmax];
    int nw, N, z[POLY_Dmax][AMBI_Dmax], m[POLY_Dmax], nz, index;
} CWS;
/*
Combined weight system: W[i][j] and d[i] are the j'th weight and the "degree"
of the i'th weight system, respectively; nw is the number of weight systems,
N is the dimension of the ambient space.
z[i][j]/m[i] are the phases of nz symmetries for polytopes on sublattices.
B describes the ambient space coordinate hyperplanes in terms of the new
(non-redundant) coordinates.
*/

typedef Long PairMat[EQUA_Nmax][VERT_Nmax];
/*
The matrix whose entries are the pairings av+c between the vertices v and
the equations (a,c).
*/

typedef struct {
    int mp, mv, np, nv, n, pic, cor, h22, h1[POLY_Dmax - 1];
} BaHo;
/*
This structure is related to Batyrev's formulas for Hodge numbers.
n     ... dimension of the polytope
pic   ... Picard number
cor   ... sum of correction terms
h1[i] ... Hodge number h_{1i}
h22   ... Hodge number h_{22} (if n = 5)
mp, mv, np, nv denote the numbers of points/vertices in the M and N lattices,
repectively.
*/

typedef struct {
    Long W[FIB_Nmax][VERT_Nmax];
    int nw, PS, ZS, nv, f[VERT_Nmax], r[VERT_Nmax], nf, nz[FIB_Nmax],
        n0[FIB_Nmax], Z[FIB_Nmax][VERT_Nmax], M[FIB_Nmax];
    GL_Long G[VERT_Nmax][POLY_Dmax][POLY_Dmax];
    PolyPointList *P;
} FibW;
/*
This list is an extension of the PolyPointList with the combined weight system.
W[i][j] is the j'th weight; nw is the number of weight systems.
*/

#ifdef __cplusplus
extern "C" {
#endif

#define INT_Nbits 32
#define LONG_LONG_Nbits 64
/*
These numbers should be set to the actual numbers of bits occupied by the
structures "unsigned int" and "unsigned long long" in your version of C.
If they are set to lower values, everything still works but may be
considerably slowed down.
*/

#if (VERT_Nmax <= INT_Nbits)
typedef unsigned int INCI;
#elif (VERT_Nmax <= LONG_LONG_Nbits)
typedef unsigned long long INCI;
#else
#define I_NUI ((VERT_Nmax - 1) / INT_Nbits + 1)
typedef struct {
    unsigned int ui[I_NUI];
} INCI;
#endif
/*
An INCI encodes the incidence relations between a face and a list of
vertices as a bit pattern (1 if a vertex lies on the face, 0 otherwise).
Depending on the allowed number VERT_Nmax of vertices, a single "unsigned int"
or "unsigned long long" may be sufficient.
If VERT_Nmax is larger than the number of bits in a "long long integer", an
array of unsigned integers is used to simulate an integer type of the required
size.
*/

typedef struct {
    int nf[POLY_Dmax + 1];              /* #(faces)[dim]  */
    INCI v[POLY_Dmax + 1][FACE_Nmax];   /*  vertex info   */
    INCI f[POLY_Dmax + 1][FACE_Nmax];   /* V-on-dual info */
    Long nip[POLY_Dmax + 1][FACE_Nmax]; /* #IPs on face  */
    Long dip[POLY_Dmax + 1][FACE_Nmax];
} FaceInfo; /* #IPs on dual  */
/*
nf[i] denotes the number of faces of dimension i
   (the number of faces of dimension n-i-1 of the dual polytope).
v[i][j] encodes the incidence relation of the j'th dim-i face with the vertices
nip[i][j] is the number of interior points of the j'th dim-i face.
f[i][j] and dip[i][j] give the same informations for the dual (n-i-1
   dimensional) faces, with f[i][j] referring to the dual vertices.
*/

#if (VERT_Nmax <= LONG_LONG_Nbits)
#define INCI_M2(x) ((unsigned char)(x) % 2)     /* value of first bit      */
#define INCI_AND(x, y) ((x) & (y))              /* bitwise logical and     */
#define INCI_OR(x, y) ((x) | (y))               /* bitwise logical or      */
#define INCI_XOR(x, y) ((x) ^ (y))              /* bitwise exclusive or    */
#define INCI_EQ(x, y) ((x) == (y))              /* check on equality       */
#define INCI_LE(x, y) INCI_EQ(INCI_OR(x, y), y) /* bitwise less or equal */
#define INCI_EQ_0(x) INCI_EQ(x, INCI_0())       /* check if all bits = 0   */
#define INCI_0() (0)                            /* set all bits to 0       */
#define INCI_1() (1)                            /* set only first bit to 1 */
#define INCI_D2(x) ((x) / 2)                    /* shift by one bit        */
#define INCI_PN(x, y) (2 * (x) + !(y))          /* shift and set first bit */
/*
For an INCI defined as a single unsigned (long long) integer whose bits are
regarded as representing incidences, these are useful definitions.
INCI_PN is particularly useful when a new vertex is added: if x represents
an equation E w.r.t. some vertex list and y is the result of evaluating E
on some new vertex V, then INCI_PN(x,y) represents x w.r.t. the vertex list
enhanced by V.
*/

#else
#define INCI_M2(x) ((x).ui[0] % 2)
INCI INCI_AND(INCI x, INCI y);
INCI INCI_OR(INCI x, INCI y);
INCI INCI_XOR(INCI x, INCI y);
int INCI_EQ(INCI x, INCI y);
int INCI_LE(INCI x, INCI y);
int INCI_EQ_0(INCI x);
INCI INCI_0();
INCI INCI_1();
INCI INCI_D2(INCI x);
INCI INCI_PN(INCI x, Long y);
#endif
/*
If we need more bits than can be represented by a single unsigned long long,
these routines are designed to simulate the above definitions.
*/

#ifdef __cplusplus
}
#endif

#endif
