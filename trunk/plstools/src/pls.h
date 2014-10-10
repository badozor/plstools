/* =============================================================================
 Title: pls.h
 Description: header file related to the file called pls.c
 Author: Pierre BADY <pierre.bady@free.fr>
 Date: 09/27/08 14:04:42
 Revision: 2010-08-27
 Version: 0.2
 Comments: RAS
 License: GPL version 2 or newer
 Copyright (C) 2008  Pierre BADY

 This program is free software; you can redistribute it and/or
 modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation; either version 2
 of the License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
============================================================================= */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <R.h>
# include <Rdefines.h>
# include <Rmath.h>
# include <Rinternals.h>
# include <R_ext/Applic.h>
# include <R_ext/Lapack.h>

// functions related to missing data
void SumX(int* n,double* X,double* S);
void ProdXY(int* n,double* x,double* y,double* S);
void SumProdXY(int* n,double* x,double* y,double* S);
void NormX(int* n,double* x);
void NormXbyY(int* n,double* x,double* y);
void replaceInteger(int n, int alpha, int* y);
void replaceDouble(int n, double alpha, double* y);
void replaceVec(int n, double* x, double* y);
void extractVec(int n, int p,int k,double* X,double* x,int mode);

// memory allocation
void taballoc(double ***tab, int l1, int c1);
void freetab(double **tab);
void vecalloc(double **vec, int n);
void freevec(double *vec);
void vecallocInt(int **vec, int n);
void freevecInt(int *vec);

// linear algebra functions
void L_prodAB(int* na, int* pa, int* nb, int* pb, double* A, double* B, double* C);
void L_prodAtB(int* na, int* pa, int* nb, int* pb, double* A, double* B, double* C);
void L_prodABt(int* na, int* pa, int* nb, int* pb, double* A, double* B, double* C);
void L_solve(int* na, int* pa, int* nb, int* pb, double* A, double* B);
void L_invA(int* na, double* A);
void L_addAB(int* na, int* pa, double* alpha, double* beta, double *A, double* B,double* C);
void L_Idn(int* na,double* A);

// plsFitter
void plsFitter(int* nc,int* nr,int* nf,double* X,double* y,double* Ch,double* Wh,double* Th,double* Ph,double* Uh);
void extractVectest(int* nr, int* nc,int* h,double* X,double* x,int* type);

// additional function used in gplsFitter
void ScaleXwt(int* n,int* p, double* X,double* wt, double* Y);
void ScaleX(int* n,int* p, double* X,double* Y);
void MeanX(int* n, double* x,double* m);
void MeanXwt(int*n, double* x,double* wt,double* m);
void varX(int* n, double* x,double* m);
void varXwt(int*n, double* x,double* wt,double* m);
void SumProdXYwt(int* n,double* x,double* y,double* wt,double* S);
void plsFitterwt(int* nc,int* nr,int* nf,double* X,double* y,double* wt, double* Ch,double* Wh,double* Th,double* Ph,double* Uh);

// additional functions realted to NIPALS functions
void rebuildTab(int n,int p, int nf, double* li, double* c1, double* A);

