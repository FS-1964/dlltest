/*******************************************************************************/
/*fftrm.c                                                                      */
/*This code contains the necessary function for fourier and inverse fourier    */
/* transformation                                                              */
/*******************************************************************************/

#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include "fftrm.h"

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

float *WR;
float *WI;

doublereal *DWR;
doublereal *DWI;

void rmpo( int *rv, int *rvp )
{
   int value_h;
   int n;

   n = 1;
   *rvp = -1;
   value_h = 1;

   while ( value_h > 0 )
   {
      value_h = *rv - n;
      (*rvp)++;
      n += n;
   }
}

void zfftrmc( doublecomplex *X, int M, int P, float D )
{
   int MV2,MM1,J,I,K,L,LE,LE1,IP,IQ,IND,IND1,R;
   int I1,J1;
   float A,B;
   float WCOS,WSIN;
   float VR,VI;
   float ARG;

   static int IPOTC;
   static float DALT;

   DWR = (doublereal *)calloc(M,sizeof(doublereal));
   DWI = (doublereal *)calloc(M,sizeof(doublereal));

   /* if (IPOTC == P & D == DALT) goto warmstart; */

   IPOTC = P;
   DALT = (float)D;
   LE = 1;
   IND = 0;

   for (L=1;L<=P;L++)
   {
      LE1 = LE;
      LE = LE*2;
      DWR[IND] = 1.0;
      DWI[IND] = 0.0;
      ARG = (float)M_PI/(float)LE1;
      WCOS = (float)cos(ARG);
      WSIN = (float)(D*sin(ARG));

      for (R=1;R<=LE1;R++)
      {
         IND1 = IND+1;
         A = (float)DWR[IND];
         B = (float)DWI[IND];
         DWR[IND1] = A*WCOS - B*WSIN;
         DWI[IND1] = B*WCOS + A*WSIN;
         ++IND;
      }
   }

   /* warmstart: */

   MV2=M/2;
   MM1=M-1;
   J=1;

   for (I=1; I<=MM1; I++)
   {
      if (I >= J)
      goto P1;

      J1 = J-1;
      I1 = I-1;

      VR = (float)RE(X[J1]);
      VI = (float)IM(X[J1]);

      RE(X[J1]) = RE(X[I1]);
      IM(X[J1]) = IM(X[I1]);

      RE(X[I1]) = VR;
      IM(X[I1]) = VI;

   P1: K = MV2;
   P2: if (K >= J) goto P3;
      J = J-K;
      K = K/2;
      goto P2;
   P3: J = J+K;
   }
   IND = 0;
   LE = 1;

   for (L=1; L<=P; L++)
   {
      LE1 = LE;
      LE = LE*2;

      for (R=0; R<LE1; R++)
      {
         WCOS = (float)DWR[IND];
         WSIN = (float)DWI[IND];
         IND = IND+1;
         for (IQ=R; IQ<M; IQ+=LE)
         {
            IP = IQ+LE1;

            A = (float)RE(X[IP]);
            B = (float)IM(X[IP]);

            VR = A*WCOS - B*WSIN;
            VI = B*WCOS + A*WSIN;

            RE(X[IP]) = RE(X[IQ]) - VR;
            IM(X[IP]) = IM(X[IQ]) - VI;

            RE(X[IQ]) = RE(X[IQ]) + VR;
            IM(X[IQ]) = IM(X[IQ]) + VI;
         }
      }
   }

   free(DWR);
   free(DWI);
}


/*=============================================================================*/
/*___1-D FFT with respect to a spatial coordinate______________________________*/
/*=============================================================================*/
int zffts( int debug, doublecomplex *X, int M )
{
   int P;
   float D;

   D = -1.0;

   rmpo( &M, &P);

   if ( debug )
   {
      printf("P = %d\n",P);
      printf("FFT ...\n");
   }

   zfftrmc( X, M, P, D); /* fftrm.c */

   return 0;
}

/*=============================================================================*/
/*___1-D Inverse FFT with respect to a spatial coordinate______________________*/
/*=============================================================================*/
int ziffts( int debug, doublecomplex *X, int M )
{
   int i;
   int P;
   float D;

   D = 1.0;

   rmpo( &M, &P);

   if ( debug )
   {
      printf("P = %d\n",P);
      printf("IFFT ...\n");
   }

   zfftrmc( X, M, P, D); /* fftrm.c */

   /*___Multiply with 1/M____*/

   for (i=0; i<M; i++)
   {
      RE(X[i]) /= (doublereal)M;
      IM(X[i]) /= (doublereal)M;
   }

   return 0;
}/*End of fftrm.c*/

