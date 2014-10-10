/* =============================================================================
 Title: pls.c
 Description: C engine for pls computation based on NIPAS algorithm
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

# include "pls.h"

/*==============================================
  fonction de Daniel Chessel pour l'allocation
  dynamique de la mémoire issue de la librairie
  ade4 (version 1.4)
==============================================*/
void taballoc (double ***tab, int l1, int c1){
    int i, j;
    if ( (*tab = (double **) calloc(l1+1, sizeof(double *))) != 0) {
        for (i=0;i<=l1;i++) {
            if ( (*(*tab+i)=(double *) calloc(c1+1, sizeof(double))) == 0 ) {
                return;
                for (j=0;j<i;j++) {
                    free(*(*tab+j));
                }
            }
        }
    }
    **(*tab) = l1;
    **(*tab+1) = c1;
}
void freetab (double **tab){
    int i, n;
    n = *(*(tab));
    for (i=0;i<=n;i++) {
            free((char *) *(tab+i) );
    }
    free((char *) tab);
}
void vecalloc (double **vec, int n){
    if ( (*vec = (double *) calloc(n+1, sizeof(double))) != 0) {
        **vec = n;
        return;
    }
    else {
        return;
    }
}
void freevec (double *vec){
    free((char *) vec);
}
void vecallocInt (int **vec, int n){
    if ( (*vec = (int *) calloc(n+1, sizeof(int))) != 0) {
        **vec = n;
        return;
    }
    else {
        return;
    }
}
void freevecInt (int *vec){
    free((char *) vec);
}
/* ============================================================
Additional functions related to linear algebra
based on libraries BLAS and LAPACK (include files in R directory)
http://www.netlib.org/lapack/
http://www.netlib.org/blas/
============================================================ */
void L_prodAB(int* na, int* pa, int* nb, int* pb, double* A, double* B, double* C){
  char *transa = "N", *transb = "N";
  double alpha = 1.0, beta = 0.0;
  F77_CALL(dgemm)(transa,transb,na,pb,pa,&alpha,A,na,B,nb,&beta,C,na);
}
void L_prodAtB(int* na, int* pa, int* nb, int* pb, double* A, double* B, double* C){
  char *transa = "T", *transb = "N";
  double alpha = 1.0, beta = 0.0;
  F77_CALL(dgemm)(transa,transb,pa,pb,na,&alpha,A,na,B,nb,&beta,C,pa);
}
void L_prodABt(int* na, int* pa, int* nb, int* pb, double* A, double* B, double* C){
  char *transa = "N", *transb = "T";
  double alpha = 1.0, beta = 0.0;
  F77_CALL(dgemm)(transa,transb,na,nb,pa,&alpha,A,na,B,nb,&beta,C,na);
}
void L_solve(int* na, int* pa, int* nb, int* pb, double* A, double* B){
  int info;
  int *ipiv;
  vecallocInt(&ipiv,*na);
  F77_CALL(dgesv)(na, pb, A, na, ipiv, B, nb, &info);
  freevecInt(ipiv);
}
void L_addAB(int* na, int* pa, double* alpha, double* beta, double *A, double* B,double* C){
  int n,p,i,j;
  n = *na;
  p = *pa;
  for(i=0;i<n;i++)
    for(j=0;j <p;j++)
      C[i+(j*n)] = (*alpha)*A[i+(j*n)] + (*beta)*B[i+(j*n)];
}
void L_Idn(int* na, double* A){
  int n,i,j;
  n = *na;
  for(i=0;i<n;i++)
    for(j=0;j <n;j++){
            if(i==j)
              A[i+(j*n)]=1.0;
            else
              A[i+(j*n)]=0.0;
    }
}
void L_invA(int* na, double* A){
  int info,n,i,j;
  int *ipiv;
  double * B;
  n= *na;
  vecallocInt(&ipiv,*na);
  vecalloc(&B,n*n);
  L_Idn(na,B);
  F77_CALL(dgesv)(na, na, A, na, ipiv, B, na, &info);
  for(i=0;i<n;i++)
    for(j=0;j <n;j++)
      A[i+(j*n)] = B[i+(j*n)];
  freevecInt(ipiv);
  freevec(B);
}
/*------------------------------------------------------------------
additional function related to the computation with missing data


------------------------------------------------------------------*/
void SumX(int* n,double* x,double* S){
  int i,nx;
  double w;
  w = 0.0;
  nx = *n;
  for(i=0;i<nx;i++)
  if(R_FINITE(x[i]))
      w +=x[i];
  S[0] = w;
}
void ProdXY(int* n,double* x,double* y,double* S){
  int i,nx;
  nx = *n;
  for(i=0;i<nx;i++)
  if(R_FINITE(x[i]) && R_FINITE(y[i]))
      S[i]=x[i]*y[i];
}
void SumProdXY(int* n,double* x,double* y,double* S){
  int i,nx;
  double w;
  w = 0.0;
  nx = *n;
  for(i=0;i<nx;i++)
  if(R_FINITE(x[i]) && R_FINITE(y[i]))
      w +=x[i]*y[i];
  S[0] = w;
}
void NormX(int* n,double* x){
  int i,nx;
  double w;
  w = 0.0;
  nx = *n;
  SumProdXY(n,x,x,&w);
  for(i=0;i<nx;i++)
    if(R_FINITE(x[i]))
      x[i] =x[i]/sqrt(w);
}
void NormXbyY(int* n,double* x,double* y){
  int i,nx;
  double w;
  w = 0.0;
  nx = *n;
  SumProdXY(n,y,y,&w);
  for(i=0;i<nx;i++)
    if(R_FINITE(x[i]))
      x[i] =x[i]/sqrt(w);
}
void replaceDouble(int n, double alpha, double* y){
  int i;
  for(i=0;i<n;i++)
    y[i] = alpha;
}
void replaceInteger(int n, int alpha, int* y){
  int i;
  for(i=0;i<n;i++)
    y[i] = alpha;
}
void replaceVec(int n, double* x, double* y){
  int i;
  for(i=0;i<n;i++)
    y[i] = x[i];
}
void extractVec(int n, int p,int k,double* X,double* x,int mode){
  int i,j;
  switch(mode){
    case 1: //extract ligne
      for(j=0;j < p;j++){
        i = k;
        x[j] = X[i+(n*j)];
      }
     break;
    default://extract columns
      for(i=0;i < n;i++){
        j = k;
        x[i] = X[i+(n*j)];
      }
    	break;
  }
}
/* plsFitter: last modified 06/14/09 15:56:49 by pbady */
void plsFitter(int* nc,int* nr,int* nf,double* X,double* y,double* Ch,double* Wh,double* Th,double* Ph,double* Uh){
  int i,h,n,p,k,m,ione;
  double auxi,S,alpha,beta,ch,dzero;
  double* yh;
  double* ytmp;
  double* th;
  double* wh;
  double* ph;
  double* uh;
  double* Xh;
  double* Xtmp1;
  double* Xtmp2;
  double* xhi;
  double* xhj;
  k = *nf;
  n = *nr;
  p = *nc;
  m = n*p;
  ch=0.0;
  alpha = 1.0;
  beta = -1.0;
  auxi = 0.0;
  S = 0.0;
  ione = 1;
  dzero = 0.0;
  vecalloc(&yh,n);
  vecalloc(&th,n);
  vecalloc(&uh,n);
  vecalloc(&wh,p);
  vecalloc(&ph,p);
  vecalloc(&Xh,m);
  vecalloc(&xhi,p);
  vecalloc(&xhj,n);
  vecalloc(&Xtmp1,m);
  vecalloc(&Xtmp2,m);
  vecalloc(&ytmp,n);
  replaceVec(m,X,Xh);
  replaceVec(n,y,yh);
  for(h=0;h<k;h++){
    ch = 0.0;
    replaceDouble(n,dzero,th);
    replaceDouble(n,dzero,ytmp);
    replaceDouble(p,dzero,ph);
    replaceDouble(p,1/sqrt(p),wh);
    replaceDouble(p,dzero,xhi);
    replaceDouble(n,dzero,xhj);
    replaceDouble(m,dzero,Xtmp1);
    replaceDouble(m,dzero,Xtmp2);
    for(i=0;i<p;i++){
      extractVec(n,p,i,Xh,xhj,0);
      auxi=0.0;
      S=0.0;
      SumProdXY(&n,yh,xhj,&auxi);
      SumProdXY(&n,yh,yh,&S);
      wh[i] = auxi/S;
    }
    NormX(&p,wh);
    for(i=0;i<n;i++){
      extractVec(n,p,i,Xh,xhi,1);
      auxi=0.0;
      S=0.0;
      SumProdXY(&p,xhi,wh,&auxi);
      SumProdXY(&p,wh,wh,&S);
      th[i] = auxi/S;
    }
    for(i=0;i<p;i++){
      extractVec(n,p,i,Xh,xhj,0);
      auxi=0.0;
      S=0.0;
      SumProdXY(&n,xhj,th,&auxi);
      SumProdXY(&n,th,th,&S);
      ph[i] = auxi/S;
    }
  L_prodABt(&n, &ione,&p,&ione,th, ph, Xtmp1);
  L_addAB(&n,&p, &alpha, &beta,Xh,Xtmp1,Xtmp2);
  replaceVec(m,Xtmp2,Xh);
  auxi = 0.0;
  S=0.0;
  SumProdXY(&n,th,yh,&auxi);
  SumProdXY(&n,th,th,&S);
  ch = auxi/S;
  for(i=0;i<n;i++)
    uh[i] = yh[i]/ch;
  auxi = (-1)*ch;
  L_addAB(&n,&ione, &alpha, &auxi,yh,th,ytmp);
// results
  Ch[h] = ch;
  for(i=0;i<p;i++){
    Ph[i+(p*h)] = ph[i];
    Wh[i+(p*h)] = wh[i];
  }
  for(i=0;i<n;i++){
    Th[i+(n*h)] = th[i];
    Uh[i+(n*h)] = uh[i];
  }
  }
  freevec(yh);
  freevec(ytmp);
  freevec(th);
  freevec(wh);
  freevec(ph);
  freevec(uh);
  freevec(Xh);
  freevec(Xtmp1);
  freevec(Xtmp2);
  freevec(xhi);
  freevec(xhj);
}
/* additional function used in gplsFitter 06/27/09 12:07:19 */
void MeanXwt(int*n, double* x,double* wt,double* m){
  int i,nx;
  double sx,sw;
  sx = 0.0;
  sw = 0.0;
  nx = *n;
  for(i=0;i<nx;i++){
    if(R_FINITE(x[i]) && R_FINITE(wt[i]))
      sx += x[i]*wt[i];
      sw += wt[i];
    }
  m[0] = sx/sw;
}
// à vérifier: faut-il utiliser MeanX ou MeanXwt?
// 2010-08-27 by pbady
void varXwt(int* n, double* x,double* wt,double* v){
  int i,nx;
  double sx,sw;
  double m[1]={0.0};
  sx = 0.0;
  sw = 0.0;
  nx = *n;
  MeanXwt(n,x,wt,m);
  for(i=0;i<nx;i++){
    if(R_FINITE(x[i]))
      sx += (x[i]-(*m))*(x[i]-(*m))*wt[i];
      sw += wt[i];
    }
  v[0] = sx/sw;
}
void MeanX(int* n, double* x,double* m){
  int i,nx;
  double sx,sw;
  sx = 0.0;
  sw = 0.0;
  nx = *n;
  for(i=0;i<nx;i++){
    if(R_FINITE(x[i]))
      sx += x[i];
      sw += 1;
    }
  m[0] = sx/sw;
}
void varX(int* n, double* x,double* v){
  int i,nx;
  double sx,sw;
  double m[1]={0.0};
  sx = 0.0;
  sw = 0.0;
  nx = *n;
  MeanX(n,x,m);
  for(i=0;i<nx;i++){
    if(R_FINITE(x[i]))
      sx += (x[i]-(*m))*(x[i]-(*m));
      sw += 1;
    }
  v[0] = sx/sw;
}
void ScaleX(int* n,int* p, double* X,double* Y){
  int i,j,col,lig;
  double*moy,*var;
  double sum,maxvar;
  sum= 0.0;
  maxvar= 0.0;
  lig=*n;
  col=*p;
  vecalloc(&moy,col);
  vecalloc(&var,col);
  for(j=0;j < col;j++){
    sum=0;
    for(i=0;i< lig;i++)
      sum+= X[i+(lig*j)];
    moy[j] = sum/(lig);
  }
  for(j=0;j < col;j++){
    sum=0;
    for(i=0;i< lig;i++)
      sum+= (X[i+(lig*j)]-moy[j])*(X[i+(lig*j)]-moy[j]);
    var[j] = sum/(lig);
  }
  for(j=0;j < col;j++)
    if(maxvar <= var[j])
      maxvar = var[j];
  for(j=0;j < col;j++)
    if(var[j] <= (0.0000001 *maxvar))
      var[j] = 1;
  for(j=0;j < col;j++)
    for(i=0;i< lig;i++)
      Y[i+(lig*j)]=(X[i+(lig*j)]-moy[j])/sqrt(var[j]);
  freevec(moy);
  freevec(var);
}
void ScaleXwt(int* n,int* p, double* X,double* wt, double* Y){
  int i,j,col,lig;
  double*moy,*var;
  double sum,sumw,maxvar;
  sum= 0.0;
  sumw= 0.0;
  maxvar=0.0;
  lig=*n;
  col=*p;
  vecalloc(&moy,col);
  vecalloc(&var,col);
  for(i=0;i< lig;i++)
      sumw+= wt[i];
  for(j=0;j < col;j++){
    sum=0;
    for(i=0;i< lig;i++)
      sum+= X[i+(lig*j)]*wt[i];
    moy[j] = sum/sumw;
  }
  for(j=0;j < col;j++){
    sum=0;
    for(i=0;i< lig;i++)
      sum+= (X[i+(lig*j)]-moy[j])*(X[i+(lig*j)]-moy[j])*wt[i];
    var[j] = sum/sumw;
  }
  for(j=0;j < col;j++)
    if(maxvar <= var[j])
      maxvar = var[j];
  for(j=0;j < col;j++)
    if(var[j] <= (0.0000001 *maxvar))
      var[j] = 1;
  for(j=0;j < col;j++)
    for(i=0;i< lig;i++)
      Y[i+(lig*j)]=(X[i+(lig*j)]-moy[j])/sqrt(var[j]);
  freevec(moy);
  freevec(var);
}
void SumProdXYwt(int* n,double* x,double* y,double* wt,double* S){
  int i,nx;
  double w;
  w = 0.0;
  nx = *n;
  for(i=0;i<nx;i++)
  if(R_FINITE(x[i]) && R_FINITE(y[i])&& R_FINITE(wt[i]))
      w +=x[i]*y[i]*wt[i];
  S[0] = w;
}
// plsFitterwt 06/27/09 12:13:36
void plsFitterwt(int* nc,int* nr,int* nf,double* X,double* y,double* wt, double* Ch,double* Wh,double* Th,double* Ph,double* Uh){
  int i,h,n,p,k,m,ione;
  double auxi,S,alpha,beta,ch,dzero;
  double* yh;
  double* ytmp;
  double* th;
  double* wh;
  double* ph;
  double* uh;
  double* Xh;
  double* Xtmp1;
  double* Xtmp2;
  double* xhi;
  double* xhj;
  k = *nf;
  n = *nr;
  p = *nc;
  m = n*p;
  ch=0.0;
  alpha = 1.0;
  beta = -1.0;
  auxi = 0.0;
  S = 0.0;
  ione = 1;
  dzero = 0.0;
  vecalloc(&yh,n);
  vecalloc(&th,n);
  vecalloc(&uh,n);
  vecalloc(&wh,p);
  vecalloc(&ph,p);
  vecalloc(&Xh,m);
  vecalloc(&xhi,p);
  vecalloc(&xhj,n);
  vecalloc(&Xtmp1,m);
  vecalloc(&Xtmp2,m);
  vecalloc(&ytmp,n);
  replaceVec(m,X,Xh);
  replaceVec(n,y,yh);
  for(h=0;h<k;h++){
    ch = 0.0;
    replaceDouble(n,dzero,th);
    replaceDouble(n,dzero,ytmp);
    replaceDouble(p,dzero,ph);
    replaceDouble(p,1/sqrt(p),wh);
    replaceDouble(p,dzero,xhi);
    replaceDouble(n,dzero,xhj);
    replaceDouble(m,dzero,Xtmp1);
    replaceDouble(m,dzero,Xtmp2);
    for(i=0;i<p;i++){
      extractVec(n,p,i,Xh,xhj,0);
      auxi=0.0;
      S=0.0;
      SumProdXYwt(&n,yh,xhj,wt,&auxi);
      SumProdXYwt(&n,yh,yh,wt,&S);
      wh[i] = auxi/S;
    }
    NormX(&p,wh);
    for(i=0;i<n;i++){
      extractVec(n,p,i,Xh,xhi,1);
      auxi=0.0;
      S=0.0;
      SumProdXY(&p,xhi,wh,&auxi);
      SumProdXY(&p,wh,wh,&S);
      th[i] = auxi/S;
    }
    for(i=0;i<p;i++){
      extractVec(n,p,i,Xh,xhj,0);
      auxi=0.0;
      S=0.0;
      SumProdXYwt(&n,xhj,th,wt,&auxi);
      SumProdXYwt(&n,th,th,wt,&S);
      ph[i] = auxi/S;
    }
  L_prodABt(&n, &ione,&p,&ione,th, ph, Xtmp1);
  L_addAB(&n,&p, &alpha, &beta,Xh,Xtmp1,Xtmp2);
  replaceVec(m,Xtmp2,Xh);
  auxi = 0.0;
  S=0.0;
  SumProdXYwt(&n,th,yh,wt,&auxi);
  SumProdXYwt(&n,th,th,wt,&S);
  ch = auxi/S;
  for(i=0;i<n;i++)
    uh[i] = yh[i]/ch;
  auxi = (-1)*ch;
  L_addAB(&n,&ione, &alpha, &auxi,yh,th,ytmp);
// results
  Ch[h] = ch;
  for(i=0;i<p;i++){
    Ph[i+(p*h)] = ph[i];
    Wh[i+(p*h)] = wh[i];
  }
  for(i=0;i<n;i++){
    Th[i+(n*h)] = th[i];
    Uh[i+(n*h)] = uh[i];
  }
  }
  freevec(yh);
  freevec(ytmp);
  freevec(th);
  freevec(wh);
  freevec(ph);
  freevec(uh);
  freevec(Xh);
  freevec(Xtmp1);
  freevec(Xtmp2);
  freevec(xhi);
  freevec(xhj);
}
//fonctions associées à la programmations de la fonction NIPALS
void rebuildTab(int n,int p, int nf, double* li, double* c1, double* A){
	int i,j,k,m;
	double dzero;
	dzero = 0.0;
	m = n*p;
	replaceDouble(m,dzero,A);
	for (k=0;k<=nf;k++)
		for(i=0;i<=n;i++)
			for(j=0;j<=n;i++)
				A[i+(j*n)] = A[i+(j*n)] + li[i+(k*nf)] * c1[j+(k*nf)];
}


