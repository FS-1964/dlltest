/************************************************************************/
/*fftrm.h                                                               */
/*This is the header file for fftrm.c                                   */
/************************************************************************/

#ifndef FFTRM_H
#define FFTRM_H

#define RE(z) ((z).r)
#define IM(z) ((z).i)

typedef float real;
typedef double doublereal;
typedef struct { real r, i; } complex;
typedef struct { doublereal r, i; } doublecomplex;

int zffts (int debug, doublecomplex *X, int M);
int ziffts (int debug, doublecomplex *X, int M);
void zfftrmc(doublecomplex *X, int M, int P, float D);
void rmpo (int *rv, int *rvp );

#endif

