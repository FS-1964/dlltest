/** * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * /
/ * * * This program extract the envelope of modulated carrier      * * * /
/ * * * Input:                                                      * * * /
/ * * * File in text format containing a table of two columns       * * * /
/ * * * (time and test PCD output voltage vd)                       * * * /
/ * * *                                                             * * * /
/ * * * Data format of input-file:                                  * * * /
/ * * *                                                             * * * /
/ * * * One data-point per line,                                    * * * /
/ * * *                                                             * * * /
/ * * * {time[seconds], sense-coil-voltage[volts])                  * * * /
/ * * *                                                             * * * /
/ * * * Data-points shall be equidistant time                       * * * /
/ * * * Minimum sampling rate: 100 MSamples/second                  * * * /
/ * * * example for spreadsheet file (start in next line):          * * * /
/ * * * (time) , (voltage )                                         * * * /
/ * * * 3.00000e-06,1.00                                            * * * /
/ * * * 3.00200e-06,1.01                                            * * * /
/ * * *                                                             * * * /
/ * * * Run:                                                        * * * /
/ * * * hilbert Filename.txt                                        * * * /
/ * * * or                                                          * * * /
/ * * * hilbert (default file name input.txt)                       * * * /
/ * * *                                                             * * * /

/*************************************************************************/
/*hilbert.c                                                              */
/*Main program                                                           */
/*************************************************************************/

# include <stdio.h>
# include <math.h>
# include <malloc.h>
# include <ctype.h>
# include <string.h>
# include "fftrm.h"

#define MAX_POINT 200000
#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

int debug=0;
int fftdebug=0;

double *Gvalue;
double *Gtime;
double *Gr;
double *Gi;
double *Gc;
doublecomplex *Gt_ifft;

/*File containing the input data*/
char *InputFileName ="input.txt";

int SampledPoints=0;
int N;
int row;
const int col=2;

int ReadData(void)
{
   float a,b;
   int i=0;
   FILE *fp1;
   i=0;
   SampledPoints=0;  //IA

   if ((fp1 = fopen(InputFileName,"r")) == NULL)
   {
      printf("Cannot open input file.\n");
      return 1;
   }

   while(!feof(fp1))
   {
      fscanf(fp1,"%e,%e\n", &a, &b);

      Gtime[SampledPoints] = a;
      Gvalue[SampledPoints] = b;
      SampledPoints++;
      if (SampledPoints>= MAX_POINT) break;
   }
   fclose(fp1);

   fp1=fopen("inputfile.txt","w");
   if (!fp1)
   {
      fprintf(stdout,"Can't write the sampled data in inputfile.txt. \n");
      return 1;
   }
   for(i=0; i<SampledPoints; i++)
   fprintf(fp1,"%e\n",Gvalue[i]); /*Gtime[i] has been omitted*/
   fclose(fp1);

   if(debug)
   {
      fp1=fopen("inputtime.txt","w");
      if (!fp1)
      {
         fprintf(stdout,"Can't write the sampled data in inputtime.txt. \n");
         return 1;
      }
      for(i=0; i<SampledPoints; i++)
      fprintf(fp1,"%e\n",Gtime[i]); /*Gtime[i] has been omitted*/
      fclose(fp1);
   }

   if(debug)
   {
      if((fp1=fopen("inputfile.bin","wb"))!=NULL)
      {
         fwrite(Gvalue,sizeof(double),SampledPoints,fp1);
         fclose(fp1);
      }
   }

   if(SampledPoints<N)
   {
      for(i=SampledPoints;i<=N;i++)
      {
         Gvalue[i] = 0;
      }
   }
   return 0;
}/*End Of Function ReadData;*/


void Fft(void)
{
   doublecomplex *Gt_freq;

   FILE *fp1,*fp2,*fp3;
   int k,num1,num2;

   Gt_freq = (doublecomplex *)calloc(sizeof(doublecomplex),row);

   /* FFT Procedure Starts for Sampled Data*/
   for(k=0;k<=N;k++)
   {
      RE(Gt_freq[k])=Gvalue[k];
      IM(Gt_freq[k])=0.0;
   }

   if(debug)
   {
      if((fp3=fopen("f.bin","wb"))!=NULL)
      {
         fwrite(Gvalue,sizeof(double),row,fp3);
         fclose(fp3);
      }
   }

   zffts(fftdebug,Gt_freq,row);/*FFT is done in spatial coordicate*/

   for (k=0;k<=N;k++)
   {
      Gr[k]=RE(Gt_freq[k]);
      Gi[k]=IM(Gt_freq[k]);
   }
   /* FFT Procedure Ends for Sampled Data*/

   /* Writing The Real And Imaginary Part Of Reflected Part for Debuging*/
   /* Writing the real part of sampled data*/

   if(debug)
   {
      if((fp1=fopen("Gr.bin","wb"))!=NULL)
      {
         num1=fwrite(Gr,sizeof(double),row,fp1);
         fclose(fp1);
      }
      else
         fprintf(stdout,"Can't Open Gr.bin");

      // Writing the img part of sampled data
      if((fp2=fopen("Gi.bin","wb"))!=NULL)
      {
         num2=fwrite(Gi,sizeof(double),row,fp2);
         fclose(fp2);
      }
      else
         fprintf(stdout,"Can't Open Gi.bin");

      fprintf(stdout,"Num of Real Part Data after FFT = %d\n",num1);
      fprintf(stdout,"Num of Img Part Data after FFT = %d\n",num2);
   }

   free(Gt_freq);

}/* End Of The Function Fft */


void PhaseShifting(void)
{
   double *tempr, *tempi;
   int k;
   FILE *fp1;

   tempr = (double *)calloc(sizeof(double),row);
   tempi = (double *)calloc(sizeof(double),row);

   for ( k=0; k<=N; k++ )
   {
      tempr[k]=Gr[k];
      tempi[k]=Gi[k];
   }

   for ( k=0; k<=ceil(N/2); k++ )
   {
      Gr[k] = tempi[k];
      Gi[k] = -tempr[k];
   }

   for ( k=(int)ceil(N/2)+1; k<=N; k++ )
   {
      Gr[k] = -tempi[k];
      Gi[k] = tempr[k];
   }

   if(debug)
   {
      if((fp1=fopen("ffrpt.bin","wb"))!=NULL)
      {
         fwrite(Gr,sizeof(double),row,fp1);
         fclose(fp1);
      }
      if((fp1=fopen("ffipt.bin","wb"))!=NULL)
      {
         fwrite(Gi,sizeof(double),row,fp1);
         fclose(fp1);
      }
   }
   free (tempr);
   free (tempi);
}/*End of PhaseShift() function*/


void Ifft(void)
{
   double *Gt_tmp; /* It takes the real part of R_ifft*/
   double *Gt_tmpi;
   FILE *fp1;
   int k,i;

   Gt_tmp = (double *)calloc(sizeof(double),row);
   Gt_tmpi = (double *)calloc(sizeof(double),row);

   for (k=0;k<=N;k++)
   {
      Gt_ifft[k].r=Gr[k];
      Gt_ifft[k].i=Gi[k];
   }

   ziffts(fftdebug,Gt_ifft,row);/*IFFT of the signal in spatial coordinate*/

   // End of IFFT
   for (k=0;k<=N;k++)
   {
      Gt_tmp[k]=Gt_ifft[k].r;
   }

   if(debug)
   {
      fp1=fopen("ifft.txt","w");
      if (!fp1)
         fprintf(stdout,"Can't write in file");
      for(i=0; i<=N; i++)
         fprintf(fp1,"%.4e\n",(Gt_ifft[i].r));
      fclose(fp1);
   }

   if(debug)
   {
      if((fp1=fopen("iffrpt.bin","wb"))!=NULL)
      {
         fwrite(Gt_tmp,sizeof(double),row,fp1);
         fclose(fp1);
      }
      if((fp1=fopen("iffipt.bin","wb"))!=NULL)
      {
         fwrite(Gt_tmpi,sizeof(double),row,fp1);
         fclose(fp1);
      }
   }
   free(Gt_tmp );
   free(Gt_tmpi );
}/* End Of Function Ifft*/


int EnvelopeReconstruction(void)
{
   FILE *fp1;
   int k;

   doublecomplex *G; /*Input signal readed from input file in complex form*/
   doublecomplex *Ganalytical;/*Analytical function of our input signal*/

   double *test;
   double *sqrtr;
   double *sqrti;

   G = (doublecomplex *)calloc(sizeof(doublecomplex),row);
   Ganalytical = (doublecomplex *)calloc(sizeof(doublecomplex),row);

   test = (double *)calloc(sizeof(double),row);
   sqrtr=(double *)calloc(sizeof(double),row);
   sqrti=(double *)calloc(sizeof(double),row);

   for (k=0;k<=N;k++)
   {
      RE(G[k]) = Gvalue[k];
      IM(G[k]) = 0.0;
   }

   for (k=0;k<=N;k++)
   {
      RE(Ganalytical[k])=G[k].r;
      IM(Ganalytical[k])=Gt_ifft[k].r;
   }

   for (k=0;k<=N;k++)
   {
      sqrtr[k]=sqrt(Ganalytical[k].r*Ganalytical[k].r+Ganalytical[k].i*Ganalytical[k].i);
   }

   fp1=fopen("output.txt","w");
   if (!fp1)
   {
      fprintf(stdout,"Can't write extracted envelope in output.txt.\n");
      free(G);
      free(Ganalytical);
      free(test);
      free(sqrtr);
      free(sqrti);
      return 1;
   }
   for(k=0; k<SampledPoints; k++)
      fprintf(fp1,"%e,%e\n",Gtime[k],sqrtr[k]);
   fclose(fp1);

   free(G);
   free(Ganalytical);
   free(test);
   free(sqrtr);
   free(sqrti);
   return 0;
}

/*Main Function*/
void hilbert(char *fnamep)
{
   int status=0,i=1;
   char fname[256];
   strcpy(fname, fnamep);
   InputFileName= fname;

   //Reading the sampled data
   do
   {
      N=(int)pow(2,i)-1;
      i++;
   }while (MAX_POINT > N);

   if (debug)
      printf("N= %d\n",N);

   row=N+1;

   Gvalue = (double *)calloc(sizeof(double),row);
   Gtime = (double *)calloc(sizeof(double),row);
   Gr = (double *)calloc(sizeof(double),row);
   Gi = (double *)calloc(sizeof(double),row);
   Gt_ifft = (doublecomplex *)calloc(sizeof(doublecomplex),row);
   Gc = (double *)calloc(sizeof(double),row);

   status = ReadData();
   if (status== 1) goto MainExit;

   /*Does FFT*/
   Fft();

   /*Appropriate Phase has been Shifted*/
   PhaseShifting();

   /*Does IFFT*/
   Ifft();

   /*Envelope Reconstruction */
   status = EnvelopeReconstruction();
   if (status== 1) goto MainExit;

   MainExit:
   free(Gvalue);
   free(Gtime);
   free(Gr);
   free(Gi);
   free(Gt_ifft);
   free(Gc);

}/*End Of Main*/

