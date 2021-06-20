#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "structures.h"
#include "hilbert.h"
#include "ovst1234.h"
#define MAX_SAMPLES 200000
#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif
/********************************************ADD new*********************************************/

/********************************************ADD new**********************************************/

void skip_line(FILE *f)
{
	int a = 0;

	a = fgetc(f);
	while (a != '\n' && a != EOF)
	{
		a = fgetc(f);
	}
}
void read_line(FILE *f, char *s, char *t)
{
	int a = 0;
	int i = 0;


	a = fgetc(f); /* Takes chars from the file pointed by f */
	while (isspace(a)) /* spaces at the beginning of line are taken out */
	{
		a = fgetc(f);
	}
	while (a != ',' && a != EOF) /* spaces at the beginning of line are taken out */
	{
		t[i++] = (char)a;
		a = fgetc(f);
	}
	t[i] = '\0';
	i = 0;
	a = fgetc(f);
	while (a != '\n' && a != EOF)
	{
		s[i++] = (char)a;
		a = fgetc(f);
	}
	s[i] = '\0'; /* We add the line end */
}

/**************************************************ADD new*********************************************************/
// Creates a ESL node with the given volt and time parameters
ESL *createnodef(double voltf, double timef)
{
	ESL *new = NULL;
	new = (ESL *)malloc(sizeof(ESL));
	if (new != NULL)
	{
		new->volt = voltf;
		new->time = timef;
		new->sig = NULL;
	}
	else
		fprintf(stderr, "Memory Error");
	return new;
}

// Frees the allocated memory for ESL nodes
void freelist(ESL *first) // frees ESL list
{
	ESL *to_free;
	while (first->sig != NULL)
	{
		to_free = first;
		first = first->sig;
		free(to_free);
	}
}

// Creates a TIME node with the given volt and time parameters
void createtime(TIMES *new, double tr, double tf, double b, double trstartind, double trendind, double tfstartind, double tfendind, double t1, double t1startind, double t1start, double t1endind, double t2, double t2startind, double t2start, double t3, double t3end, double t3endind, double t4, double t4endind, double t5, double t5startind, double t6, double t6end, double t6endind, double a, double tploone)
{
	new->tf = tf;
	new->tr = tr;
	new->b = b;
	new->trstartind = trstartind;
	new->trendind = trendind;
	new->tfstartind = tfstartind;
	new->tfendind = tfendind;
	new->t1 = t1;
	new->t1startind = t1startind;
	new->t1start = t1start;
	new->t1endind = t1endind;
	new->t2 = t2;
	new->t2startind = t2startind;
	new->t2start = t2start;
	new->t3 = t3;
	new->t3end = t3end;
	new->t3endind = t3endind;
	new->t4 = t4;
	new->t4endind = t4endind;
	new->t5 = t5;
	new->t5startind = t5startind;
	new->t6 = t6;
	new->t6end = t6end;
	new->t6endind = t6endind;
	new->a = a;
	new->tploone = tploone;
}

// Inserts a ESL node (new) in a list pointed by "first"
void insert_node(ESL **first, ESL *new)
{
	ESL *p = NULL;
	ESL *previous = NULL;
	if (new != NULL)
	{
		p = *first;
		if (p == NULL)
		{
			*first = new;
		}
		else
		{
			while (p != NULL)
			{
				previous = p;
				p = p->sig;
			}
			previous->sig = new;
		}
	}
}

/* Multiplies order polynomials supposing (x^2 + b*x + c) */
/* b and c are complex values stored in a table where even elements are */
/* real and odd elements imaginary */
double *mult_poli(int num_pol, double *b, double *c)
{
	int i = 0;
	int y = 0;
	double *real;
	double *imag;
	double *vector;
	double *new_real;
	double *new_imag;
	double real_b = 0;
	double real_c = 0;
	double imag_b = 0;
	double imag_c = 0;

	real = (double *)calloc(4 * num_pol, sizeof(double));
	imag = (double *)calloc(4 * num_pol, sizeof(double));
	new_real = (double *)calloc(4 * num_pol, sizeof(double));
	new_imag = (double *)calloc(4 * num_pol, sizeof(double));
	vector = (double *)calloc(4 * num_pol, sizeof(double));

	real[0] = c[0];
	real[1] = b[0];
	real[2] = 1;
	imag[0] = c[1];
	imag[1] = b[1];
	imag[2] = 0;

	for (i = 3; i<(4 * num_pol); i++)
	{
		real[i] = 0;
		imag[i] = 0;
	}

	for (y = 1; y<num_pol; y++)
	{
		// Selects values b and c
		real_b = b[2 * y];
		real_c = c[2 * y];
		imag_b = b[2 * y + 1];
		imag_c = c[2 * y + 1];

		for (i = 0; i <= (2 * num_pol - 2); i++)
		{
			// Starts with coeff "c"
			new_real[i] += real[i] * real_c - imag[i] * imag_c;
			new_imag[i] += real[i] * imag_c + imag[i] * real_c;
			// Continues with coeff "b"
			new_real[i + 1] += real[i] * real_b - imag[i] * imag_b;
			new_imag[i + 1] += real[i] * imag_b + imag[i] * real_b;
			// Finishes with coeff "1"
			new_real[i + 2] += real[i];
			new_imag[i + 2] += imag[i];
		}

		// Update Values
		for (i = 0; i<(4 * num_pol); i++)
		{
			real[i] = new_real[i];
			imag[i] = new_imag[i];
			new_real[i] = 0;
			new_imag[i] = 0;
		}
	}

	for (y = 0; y<(2 * num_pol); y++)
	{
		vector[2 * y] = imag[y];
		vector[2 * y + 1] = real[y];
	}

	free(new_imag);
	free(new_real);
	free(imag);
	free(real);

	return (vector);
}

// Part of the calculation of the butterworth coeffs.
double *butter_d_coeffs(double freq1, double freq2)
{
	int butter_order = 2;
	int index = 0;
	double theta = 0;     // M_PI *(freq2-freq1)/2.0
	double cp = 0;        // cosine of phi
	double *vec_r = 0;    // z^-2 coefficients
	double *vec_t = 0;    // z^-1 coefficients
	double *dcoeff = 0;   // d coefficients
	double pole_ang = 0;      // pole angle
	double divisor = 0;

	cp = cos(M_PI*(freq2 + freq1) / 2.0);
	theta = M_PI*(freq2 - freq1) / 2.0;

	vec_r = (double *)calloc(2 * butter_order, sizeof(double));
	vec_t = (double *)calloc(2 * butter_order, sizeof(double));

	for (index = 0; index<butter_order; ++index)
	{
		pole_ang = M_PI*(double)(2 * index + 1) / (double)(2 * butter_order);
		divisor = sin(2 * theta)*sin(pole_ang) + 1.0;
		vec_r[2 * index] = cos(2 * theta) / divisor;
		vec_r[2 * index + 1] = sin(2 * theta)*cos(pole_ang) / divisor;
		vec_t[2 * index] = -2.0*cp*(cos(theta) + sin(theta)*sin(pole_ang)) / divisor;
		vec_t[2 * index + 1] = -2.0*cp*sin(theta)*cos(pole_ang) / divisor;
	}

	dcoeff = mult_poli(butter_order, vec_t, vec_r);

	dcoeff[4] = dcoeff[1];
	dcoeff[3] = dcoeff[3];
	dcoeff[2] = dcoeff[5];
	dcoeff[1] = dcoeff[7];
	dcoeff[0] = 1;

	for (index = 5; index <= 2 * butter_order; index++)
		dcoeff[index] = 0;

	free(vec_t);
	free(vec_r);

	return(dcoeff);
}

// Calculates the Butterworth filter coefficients
void butterworth_coeffs(double freq1, double freq2, double *dfiltercoeff, double *cfiltercoeff)
{
	// n filter order
	// freq1, freq2 lower/uppercutoff frequencies

	double sf;        	// scaling factor
	double *dcoeff;   	// d coefficients
	double cotan = 0;     // cotangent of theta

	/* calculate the d coefficients */
	dcoeff = butter_d_coeffs(freq1, freq2);

	/* d coefficients for 4th order butterworth */
	dfiltercoeff[0] = dcoeff[0];  // Always 1
	dfiltercoeff[1] = dcoeff[1];
	dfiltercoeff[2] = dcoeff[2];
	dfiltercoeff[3] = dcoeff[3];
	dfiltercoeff[4] = dcoeff[4];

	/* scalling factor for the c filter coefficients (Butterworth 4th order */
	cotan = 1.0 / tan(M_PI*(freq2 - freq1) / 2.0);
	sf = (1.0 / (((cotan + sqrt(2) / 2)*(cotan + sqrt(2) / 2)) + 1 / 2));

	/* c coefficients for 4th order butterworth*/
	cfiltercoeff[0] = 1 * sf;
	cfiltercoeff[1] = 0 * sf;
	cfiltercoeff[2] = -2 * sf;
	cfiltercoeff[3] = 0 * sf;
	cfiltercoeff[4] = 1 * sf;

	free(dcoeff);
}

// Checks if the data input is adequate to our algorithms
int datacheck(int posval, int negval, int samplesp, double tlast, FILE *pointfile)
{
	double diffr = 0.0;
	char timestr1[25000];
	char timestr2[25000];
	char voltstr[25000];
	double timestrf1 = 0;
	double timestrf2 = 0;
	double cut_sample = 0;
	double delta_t = 0;
	double val_t = 0;
	int loop = 0;
	double linf = 0;
	int lind = 0;

	// Checks that there are (nearly) much positive as negative values
	if (posval>negval)
		diffr = (posval - negval) / ((posval + negval) / 2);
	else
		diffr = (negval - posval) / ((posval + negval) / 2);

	if (diffr>0.8)
		fprintf(stdout, "Data Corrupted: Too little negative (or positive) values\n");

	// L=n*p with P=2*pi and n=1,2,3... - Cuts data
	rewind(pointfile);
	read_line(pointfile, voltstr, timestr1);
	read_line(pointfile, voltstr, timestr1);
	read_line(pointfile, voltstr, timestr1);      // Skips csv header if present or not

	read_line(pointfile, voltstr, timestr1);
	while (voltstr[0] != '\0')
	{
		read_line(pointfile, voltstr, timestr2);
		loop++;
	}
	loop = loop + 3;

	rewind(pointfile);
	for (lind = 0; lind<loop; lind++)
	{
		read_line(pointfile, voltstr, timestr2);
	}

	timestrf1 = atof(timestr1);         // t(4)
	timestrf2 = atof(timestr2);         // t(end)
	cut_sample = ((1 / 13.56e6) / ((timestrf2 - timestrf1) / (loop - 1)));
	linf = samplesp;
	while (linf>0)
	{
		linf = linf - cut_sample;
	}
	linf += cut_sample;
	samplesp = samplesp - linf - 3;

	// At least 7 points per sample
	delta_t = tlast - timestrf1;
	val_t = delta_t / samplesp;
	if (val_t>(1 / 13.56e6) / 7)
		fprintf(stdout, "More samples points needed - Nyquist\n");

	return (samplesp);
}

// Finds the most frequent value(s) of the given signal, Hmax (types A/B) and Hmin (type B).
void Hmaxfinder(double *env, double *Hmax, double *Hmin, int numsamples)
{
	int hist[2001] = { 0 }; // IA Changed memory access violation. Increased +1
	int hi_low_i = 0;
	double tophist = 0;
	double bothist = 10;
	double diffhist = 0;
	double value = 0;
	int histind = 0;
	int max_i = 0;
	int min_i = 0;
	double max = 0;
	double min = 0;

	// Finds higher and lower values of samples
	for (hi_low_i = 0; hi_low_i<MAX_SAMPLES; hi_low_i++)
	{
		if (env[hi_low_i] != 0)
		{
			if (env[hi_low_i]<bothist)
				bothist = env[hi_low_i];
			if (env[hi_low_i]>tophist)
				tophist = env[hi_low_i];      // Finds limits for the histogram
		}
	}
	diffhist = tophist - bothist;

	for (hi_low_i = 0; hi_low_i<numsamples; hi_low_i++)
	{
		if (env[hi_low_i] != 0)
		{
			value = env[hi_low_i];
			histind = (int)(2000 * ((value - bothist) / diffhist));   // Performs a lineal quantization
			hist[histind]++;
		}
	}

	for (hi_low_i = 0; hi_low_i<1000; hi_low_i++)
	{
		if (hist[hi_low_i]>min)    // Searchs most frequent value in the lower half of the form
		{
			min = hist[hi_low_i];
			min_i = hi_low_i;
		}
		*Hmin = (bothist + (diffhist / 2000)*(min_i));
	}

	for (hi_low_i = 1001; hi_low_i<2000; hi_low_i++)
	{
		if (hist[hi_low_i]>max)    // Search most frequent value in the upper half of the form
		{
			max = hist[hi_low_i];
			max_i = hi_low_i;
		}
		*Hmax = (bothist + (diffhist / 2000)*(max_i));
	}
}

// Linear convolution (z= x convolve y)
void LinearConvolution(double X[], double Y[], double Z[], int lenx, int leny)
{
	double *zptr, s, *xp, *yp;
	int lenz;
	int i, n, n_lo, n_hi;

	lenz = lenx + leny - 1;
	zptr = Z;

	for (i = 0; i<lenz; i++)
	{
		s = 0.0;
		n_lo = 0>(i - leny + 1) ? 0 : i - leny + 1;
		n_hi = lenx - 1<i ? lenx - 1 : i;
		xp = X + n_lo;
		yp = Y + i - n_lo;

		for (n = n_lo; n <= n_hi; n++)
		{
			s += *xp * *yp;
			xp++;
			yp--;
		}

		*zptr = s;
		zptr++;
	}
}

int envfilt(double *output, double *toutput, int filterlength, double tini, double tend, int lengthp, double *envelope)
{
	int cofpi = 0;
	int xx = 0;
	double cofp = 0.0;
	int lengthp1 = 0;
	double lengthf = 0;
	double cof[2000] = { 0 };
	double points = 0.0;
	int pointsi = 0;
	int lengthtotal = 0;

	cofp = (73.75e-9) / ((tend - tini) / (lengthp));
	cofpi = cofp + 0.5;
	lengthf = cofpi*filterlength;
	points = (5 * 73.75e-9) / ((tend - tini) / (lengthp - 1));
	pointsi = (int)points + 1;
	lengthp1 = lengthp;

	for (xx = 0; xx<lengthf; xx++)
		cof[xx] = 1 / lengthf;

	for (xx = lengthf + 1; xx<2000; xx++)
		cof[xx] = 0;

	LinearConvolution(cof, output, envelope, lengthf, lengthp);

	for (xx = 0; xx<(pointsi); xx++) // "Cuts" envelope
	{
		envelope[xx] = 0.0;
		toutput[xx] = 0.0;
		envelope[lengthp1 - xx] = 0.0;
		toutput[lengthp1 - xx] = 0.0;
	}

	for (xx = lengthp1 + 1; xx<MAX_SAMPLES; xx++)
	{
		envelope[xx] = 0.0;
		toutput[xx] = 0.0;
	}

	lengthtotal = lengthp1 - 2 * (pointsi);
	return (lengthtotal);
}

// Performs the search of a certain level (target) in the envelope, i.e. 5% ,60%, 90% in type A, 106 Kbps
int localizador(double *env, double *toutput, double target, ESL **crosses, int env_length)
{
	int flag = 0;
	double diff;
	ESL *new;
	double v;
	double t;
	double lasttime;
	double lastvolt;
	int crosscounter = 0;
	int locat_index = 0;
	int locat_index_start = 0;

	while (env[locat_index] == 0.0)  // Leaves 0s out
		locat_index++;

	locat_index_start = locat_index;
	if (env[locat_index] - target>0)
	{
		for (locat_index = locat_index_start; locat_index<env_length + locat_index_start - 1; locat_index++)  // IFX
		{
			diff = env[locat_index] - target;
			if (diff<0 && flag == 0 && env[locat_index] != 0.0)  // At the beginning or after an odd occurrence, envelope is over "target" level
			{
				flag = 1; // down!
				v = env[locat_index];
				t = toutput[locat_index];
				new = createnodef(v, t);
				insert_node(crosses, new);
				crosscounter++;
			}
			if (diff>0 && flag == 1 && env[locat_index] != 0.0)   // After first (or even) occurrence, envelope is under "target" level
			{
				flag = 0; // up!
				new = createnodef(lastvolt, lasttime);
				insert_node(crosses, new);
				crosscounter++;
			}
			lasttime = toutput[locat_index];
			lastvolt = env[locat_index];
		}                                       // Returns all occurrences with time and volt level in a list
	}
	else
		fprintf(stdout, "Signal is not ---|___|--- \n");

	return (crosscounter);                           // Also returns how many occurrences appeared
}

// Function that calculates the relevant times
void tfinder(char type, double *env, double *toutput, double tini, double Hmax, double Hmin, int rate, int env_length, TIMES *timeres)
{
	double *envc = NULL;
	double *envc2 = NULL;
	ESL *crosses = NULL;
	ESL *crosses2 = NULL;
	ESL *crosses3 = NULL;
	ESL *crosses_WORK = NULL;
	double ninety = 0.0;
	double five = 0.0;
	double sixty = 0.0;
	double tp90one = 0.0;
	double tp90two = 0.0;
	double tp5one = 0.0;
	double tp5two = 0.0;
	double tp60two = 0.0;
	double vp90one = 0.0;
	double vp90two = 0.0;
	double vp5one = 0.0;
	double vp5two = 0.0;
	double vp60two = 0.0;
	double tphione = 0.0;
	double tphitwo = 0.0;
	double vphione = 0.0;
	double vphitwo = 0.0;
	double tpmidone = 0.0;
	double vpmidone = 0.0;
	double tploone = 0.0;
	double tplotwo = 0.0;
	double vplotwo = 0.0;
	double t1 = 0.0;
	double t2 = 0.0;
	double t3 = 0.0;
	double t4 = 0.0;
	double t5 = 0.0;
	double t6 = 0.0;
	int flag = 0;
	int flag2 = 0;
	int flag3 = 0;
	int flag_improv = 0;
	int x_improv = 0;
	double minvolt = 0.0;
	double highrate_low = 0.0;
	double highrate_mid = 0.0;
	double highrate_hi = 0.0;
	double a = 0.0;
	double t6end = 0.0;
	double t5startind = 0.0;
	double t6endind = 0.0;
	double b = 0.0;
	double B_low = 0.0;
	double B_hi = 0.0;
	double tr = 0.0;
	double tf = 0.0;
	double tfstartind = 0.0;
	double tfendind = 0.0;
	double trstartind = 0.0;
	double trendind = 0.0;
	double t2startind = 0.0;
	double t2start = 0.0;
	double t1startind = 0.0;
	double t1start = 0.0;
	double t3end = 0.0;
	double t4endind = 0.0;
	double t3endind = 0.0;
	double t1endind = 0.0;
	double oscmin = 0.0;
	double osctmin = 0.0;
	double oscmax = 0.0;
	double osctmax = 0.0;
	ESL *crossescopy = NULL;
	double tim3 = 0.0;
	int index_A = 0;
	int index_A2 = 0;
	int index_chain = 0;
	int i = 0;

	b = Hmin;
	envc = env;
	envc2 = envc;

	switch (type)
	{
	case 'A':
	{
		switch (rate)
		{
		case 106:
		{
			ninety = Hmax*0.9;
			five = Hmax*0.05;
			sixty = Hmax*0.6;

			flag2 = localizador(envc, toutput, five, &crosses2, env_length);      // Finds 5% of Hmax
			if (flag2 == 2)      // if there are two occurrences, there큦 no problem...
			{
				tp5one = crosses2->time;         // Temporary values are stored for future use
				vp5one = crosses2->volt;
				tp5two = crosses2->sig->time;
				vp5two = crosses2->sig->volt;
				freelist(crosses2);
			}
			else if (flag2 == 0)          // ...if there is no occurrence...
				fprintf(stdout, "5 percent of Hmax not reached - maybe wrong type or bitrate? \n");

			else if (flag2>2)         // ...if there are more than two occurrences...
			{                     // ...it must be checked that "peaks" comply the ISO restrictions

				while (toutput[index_A]<crosses2->sig->time)
					index_A++;

				oscmin = envc2[index_A];
				for (index_chain = 0; index_chain<flag2 - 2; index_chain++)
				{
					crosses2 = crosses2->sig;
				}

				while (toutput[index_A] <= crosses2->time)
				{
					if (envc2[index_A]>oscmax)
					{
						oscmax = envc2[index_A];
						osctmax = toutput[index_A];
					}
					index_A++;
				}

				while (envc2[index_A2] == 0)
				{
					index_A2++;
				}

				while (envc2[index_A2]>oscmax)
				{
					index_A2++;
				}

				oscmin = envc2[index_A2];
				osctmin = toutput[index_A2];

				if (osctmax - osctmin>5e-7)
					fprintf(stdout, "Monotony not fulfilled \n");

				tp5one = crosses2->time;         // Temporary values are stored for future use
				vp5one = crosses2->volt;
				tp5two = crosses2->sig->time;
				vp5two = crosses2->sig->volt;
				freelist(crosses2);
			}

			flag = localizador(envc, toutput, ninety, &crosses, env_length);      // Finds 90% of Hmax
			if (flag >= 2)
			{
				crosses_WORK = crosses;   // Copy of crosses to work with
				while (x_improv<flag)
				{
					if (crosses_WORK->time<tp5one)
					{
						tp90one = crosses_WORK->time;      // Temporary values are stored for future use
						vp90one = crosses_WORK->volt;
					}

					if (crosses_WORK->time>tp5two && flag_improv == 0)
					{
						tp90two = crosses_WORK->time;      // Temporary values are stored for future use
						vp90two = crosses_WORK->volt;
						flag_improv = 1;
					}
					crosses_WORK = crosses_WORK->sig;
					x_improv++;
				}
			}

			else             // ...otherwise...
			{
				fprintf(stdout, "90 %% of Hmax not found - Noise Too High \n");
			}

			freelist(crosses);

			flag3 = localizador(envc, toutput, sixty, &crosses3, env_length); // Finds 60% of Hmax
			if (flag3 == 2)               // if there are two occurrences, there큦 no problem...
			{
				tp60two = crosses3->sig->time;      // Temporary values are stored for future use
				vp60two = crosses3->sig->volt;
				freelist(crosses3);
			}

			t1 = tp5two - tp90one;      // Definitive values are calculated and stored for display
			t2 = tp5two - tp5one;
			t3 = tp90two - tp5two;
			t4 = tp60two - tp5two;

			t1start = tp90one;      // Other important values for the coming functions
			t2start = tp5one;
			t3end = tp90two;
			t1startind = vp90one;
			t1endind = vp5two;
			t2startind = vp5one;
			t3endind = vp90two;
			t4endind = vp60two;

			createtime(timeres, 0, 0, 0, 0, 0, 0, 0, t1, t1startind, t1start, t1endind, t2, t2startind, t2start, t3, t3end, t3endind, t4, t4endind, 0, 0, 0, 0, 0, 0, 0);

		}
		break;

		case 212:
		case 424:
		case 848:
		{
			ninety = Hmax*0.9;
			while (env[index_A] == 0.0)        // Finds first value different of 0.0
				index_A++;

			minvolt = env[index_A];
			while (env[index_A] != 0.0)        // All values are considered
			{
				if (env[index_A]<minvolt)
				{
					minvolt = env[index_A];   // Finds minimal voltage
				}
				index_A++;
			}

			highrate_low = minvolt + 0.1*(Hmax - minvolt);   // Calculates target
			flag = localizador(envc, toutput, highrate_low, &crosses, env_length);  // Finds target
			if (flag == 2)               // ...if there are two occurrences, there큦 no problem...
			{
				tploone = crosses->time;         // Temporary values are stored for future use
				tplotwo = crosses->sig->time;
				vplotwo = crosses->sig->volt;
			}

			else if (flag>2)         // if there are more than two occurrences...
			{                     // ...it must be checked that "peaks" comply the ISO restrictions
				while (toutput[index_A2]<crosses->time || toutput[index_A2] == 0)
					index_A2++;

				//if envc2[indexA2]!=0;
				oscmin = envc2[index_A2];
				while (envc2[index_A2] <= oscmin)
				{
					oscmin = envc2[index_A2];
					index_A2++;
				}
				osctmin = toutput[index_A2];

				crossescopy = crosses;
				for (i = 1; i<(flag - 1); i++)
					crossescopy = crossescopy->sig;

				tim3 = crossescopy->time;
				while (toutput[index_A2]<tim3)
					index_A2++;

				oscmax = envc2[index_A2];
				while (toutput[index_A2]<tim3)
				{
					if (oscmax<envc2[index_A2])
					{
						oscmax = envc2[index_A2];
						osctmax = toutput[index_A2];
					}
					index_A2++;
				}

				if (oscmax - oscmin>(0.09*(Hmax - oscmin)))
					fprintf(stdout, "Monotony not fulfilled \n");

				for (i = 1; i<(flag - 1); i++)
					crosses = crosses->sig;

				tploone = crosses->time;         // Temporary values are stored for future use
				tplotwo = crosses->sig->time;
				vplotwo = crosses->sig->volt;
			}
			freelist(crosses);

			highrate_hi = ninety + 0.1*minvolt;      // Calculates target
			flag = localizador(envc, toutput, highrate_hi, &crosses2, env_length);   // Finds target
			if (flag >= 2)
			{
				crosses_WORK = crosses2;
				while (x_improv<flag)
				{
					if (crosses_WORK->time<tploone)
					{
						tphione = crosses_WORK->time;      // Temporary values are stored for future use
						vphione = crosses_WORK->volt;
					}

					if (crosses_WORK->time>tplotwo && flag_improv == 0)
					{
						tphitwo = crosses_WORK->time;      // Temporary values are stored for future use
						vphitwo = crosses_WORK->volt;
						flag_improv = 1;
					}
					crosses_WORK = crosses_WORK->sig;
					x_improv++;
				}
			}

			else                         // ...otherwise...
			{
				fprintf(stdout, "90 %% of Hmax not reached! - Noise Too High? \n");
			}

			freelist(crosses2);

			highrate_mid = (Hmax + minvolt) / 2;      // Calculates target
			flag = localizador(envc, toutput, highrate_mid, &crosses3, env_length);   // Finds target
			if (flag == 2)                  // ...if there are two occurrences, there큦 no problem...
			{
				tpmidone = crosses3->time;      // Temporary values are stored for future use
				vpmidone = crosses3->volt;
				freelist(crosses3);
			}
			else                         // ...otherwise...
				fprintf(stdout, "Noise Too High \n");

			t1 = tplotwo - tphione;         // Definitive values are calculated and stored for display
			t5 = tplotwo - tpmidone;
			t6 = tphitwo - tplotwo;
			a = minvolt;

			t6end = tphitwo;            // Other important values for the coming functions
			t1start = tphione;
			t1startind = vphione;
			t5startind = vpmidone;
			t1endind = vplotwo;
			t6endind = vphitwo;

			createtime(timeres, 0, 0, 0, 0, 0, 0, 0, t1, t1startind, t1start, t1endind, 0, 0, 0, 0, 0, 0, 0, 0, t5, t5startind, t6, t6end, t6endind, a, tploone);
		}
		break;
		}
	}
	break;

	case 'B':
	{
		B_low = b + 0.1*(Hmax - b);         // Calculates target
		flag = localizador(envc, toutput, B_low, &crosses, env_length);   // Finds target
		if (flag >= 2)
		{
			crosses_WORK = crosses;
			tploone = crosses_WORK->time;            // Temporary values are stored for future use
			while (x_improv<flag)
			{
				tplotwo = crosses_WORK->time;         // Temporary values are stored for future use
				vplotwo = crosses_WORK->volt;
				crosses_WORK = crosses_WORK->sig;
				x_improv++;
			}
		}
		else
		{
			fprintf(stdout, "Monotony not fulfilled\n");
		}

		freelist(crosses);

		B_hi = Hmax - 0.1*(Hmax - b);         // Calculates target
		flag = localizador(envc, toutput, B_hi, &crosses2, env_length);   // Finds target
		if (flag >= 2)
		{
			x_improv = 0;
			flag_improv = 0;
			crosses_WORK = crosses2;
			while (x_improv<flag)
			{
				if (crosses_WORK->time<tploone)
				{
					tphione = crosses_WORK->time;      // Temporary values are stored for future use
					vphione = crosses_WORK->volt;
				}

				if (crosses_WORK->time>tplotwo && flag_improv == 0)
				{
					tphitwo = crosses_WORK->time;      // Temporary values are stored for future use
					vphitwo = crosses_WORK->volt;
					flag_improv = 1;
				}
				crosses_WORK = crosses_WORK->sig;
				x_improv++;
			}
		}

		else
		{
			fprintf(stdout, "Monotony not fulfilled\n");
		}

		freelist(crosses2);

		tf = tploone - tphione;      // Definitive values are calculated and stored for display
		tr = tphitwo - tplotwo;
		tfstartind = tphione;      // Other important values for the coming functions
		tfendind = tploone;
		trstartind = tplotwo;
		trendind = tphitwo;

		createtime(timeres, tr, tf, b, trstartind, trendind, tfstartind, tfendind, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
	}
	break;
	}
}

// Checks monotony on the falling edge
void monocheck(double *env, double *toutput, double Hmax, TIMES *timesp, int rate, char type)
{
	double tinit = 0.0;
	double tend = 0.0;
	double compare = 0.0;
	double timer0 = 0.0;
	double timer1 = 0.0;
	double volt0 = 0.0;
	double volt1 = 0.0;
	int counter = 0;
	int flag_mono = 0;

	switch (type)
	{
	case 'A':
	{
		switch (rate)
		{
		case 106:
		{
			while (env[counter] == 0)
				counter++;

			tinit = timesp->t1start;
			tend = timesp->t2start;

			while (toutput[counter]<tinit)
				counter++;  //find first value

			while (toutput[counter]<tend)
			{
				compare = env[counter];
				if (compare<env[counter + 1])
				{
					timer0 = toutput[counter];
					volt0 = env[counter];
					while (volt0<env[counter + 1]) // growing values...
					{
						counter++;
						volt0 = env[counter];
					}
					timer1 = toutput[counter]; // ...max.value -> time
					if (timer1 - timer0>5e-6)
						fprintf(stdout, "Monotony not fulfilled \n");
				}
				else
					counter++;
			}
		}
		break;

		case 212:
		case 424:
		case 848:
		{
			while (env[counter] == 0)
				counter++;

			tinit = timesp->t1start;
			tend = timesp->t1endind;
			while (toutput[counter]<tinit)
				counter++;  //find first value

			while (env[counter]>tend)
			{
				compare = env[counter];
				if (compare<env[counter + 1])
				{
					volt0 = env[counter];
					volt1 = volt0;
					while (volt0<env[counter + 1]) // growing values...
					{
						counter++;
						volt0 = env[counter];
					}
					if (volt1 - volt0>0.09*(Hmax - volt0))
						fprintf(stdout, "Monotony not fulfilled \n");
				}
				else
					counter++;
			}
		}
		break;
		}
	}
	break;

	case 'B':
	{
		while (env[counter] == 0)
			counter++;

		tinit = timesp->tfstartind;
		tend = timesp->tfendind;
		while (toutput[counter]<tinit)
			counter++;  //find first value

		while (toutput[counter]<tend)
		{
			compare = env[counter];
			if (compare<env[counter++])
				flag_mono = 1;
		}
		if (flag_mono == 1)
			fprintf(stdout, "Monotony not fulfilled \n");
	}
	break;
	}
}

// Function that calculates the overshoot times
void overshoot(TIMES *timesp, double Hmax, double *env2, double *toutput, int rate, char type, int samples, SHOOTREADER *shootreader)
{
	double shootind = 0.0;
	double shootind_b = 0.0;
	double hr_reader = 0.0;
	double hf_reader = 0.0;
	double above = 0.0;
	double above_b = 10.0;
	double start = 0.0;
	int index_samples = 0;

	switch (type)
	{
	case 'A':
	{
		switch (rate)
		{
		case (106) :
		{
			start = timesp->t3end;
			while (toutput[index_samples] <= start)
				index_samples++;

			while (index_samples <= samples)
			{
				if (env2[index_samples]>above)
				{
					above = env2[index_samples];
					shootind = toutput[index_samples];
				}
				index_samples++;
			}
		}
				   break;
		case (212) :
		case (424) :
		case (848) :
		{
			start = timesp->t6end;
			while (toutput[index_samples] <= start)
				index_samples++;

			while (index_samples <= samples)
			{
				if (env2[index_samples]>above)
				{
					above = env2[index_samples];
					shootind = toutput[index_samples];
				}
				index_samples++;
			}

			if (above<Hmax)      // In very strange cases if there큦 no overshoot, the highest point
				above = Hmax;      // in the curve can be cutted off by envfilt, producing a negative hr

		}
				   break;
		}
	}
	break;
	case 'B':
	{
		start = timesp->trendind;
		while (toutput[index_samples]<start)         // Starts at the rising edge
			index_samples++;

		while (index_samples <= samples)
		{
			if (env2[index_samples]>above)
			{
				above = env2[index_samples];
				shootind = toutput[index_samples];
				hr_reader = (above - Hmax) / (Hmax - timesp->b);
				if (hr_reader<0)   // In very strange cases if there큦 no overshoot, the highest point
					hr_reader = 0;   // in the curve can be cutted off by envfilt, producing a negative hr
			}
			index_samples++;
		}

		index_samples = 0;
		start = timesp->tfendind;
		while (toutput[index_samples] == 0)
			index_samples++;
		while (toutput[index_samples]<start)
			index_samples++;

		while (toutput[index_samples]<(timesp->trstartind))
		{
			if (env2[index_samples]<above_b && env2[index_samples] != 0)
			{
				above_b = env2[index_samples];
				shootind_b = toutput[index_samples];
				hf_reader = (timesp->b - above_b) / (Hmax - timesp->b);
			}
			index_samples++;
		}
	}
	break;

	}
	shootreader->shootind = shootind;
	shootreader->shootind_b = shootind_b;
	shootreader->hr_reader = hr_reader;
	shootreader->hf_reader = hf_reader;
	shootreader->above = above;
	shootreader->above_b = above_b;
}

// Calculates the modulation index "m"
double modulation(char type, double Hmax, double b)
{
	double m = 0;
	switch (type)
	{
	case 'A':
	{
		// m is not defined for type A
	}
	break;
	case 'B':
	{
		m = 100 * (Hmax - b) / (Hmax + b); // In %
	}
	break;
	}
	return (m);
}

// Displays on screen the results of the calculations
OVST1234 display(char type, int rate, SHOOTREADER *shootreader2, TIMES *timesp, double Hmax, double m, double* myovst1234)
{
	double ovs = 0;
	double ovsb1 = 0;
	double ovsb2 = 0;
	OVST1234 my_ovst1234;
	my_ovst1234.ovs = 0.0;
	my_ovst1234.t1 = 0.0;
	my_ovst1234.t2 = 0.0;
	my_ovst1234.t3 = 0.0;
	my_ovst1234.t4 = 0.0;

	//  fprintf (stdout,"\n");    // 2nd set of functions, on debug purposes
	switch (type)
	{
	case 'A':
	{

		/* fprintf(stdout,"---RESULTS------------\n");
		fprintf(stdout,"Type A - Bitrate %d\n", rate);
		fprintf(stdout,"---Overshoot----------\n");*/
		ovs = (((shootreader2->above) - Hmax) / (Hmax - timesp->a)) * 100;			//<== OVS should be displayed

		myovst1234[0] = ovs; /*** changeart= Added save the ovs value to export ***/
		my_ovst1234.ovs = ovs;
		//if (ovs>0)
		//   fprintf(stdout,"Overshoot = %f %% \n", ovs);     /*** changeart= commented ***/
		//else
		//   fprintf(stdout,"Overshoot = 0 %% \n");         /*** changeart= commented ***/
		switch (rate)
		{
		case (106) :
		{
			/* fprintf(stdout,"---timings------------\n");*/
			/* fprintf(stdout,"t1 = %f microsec. \n", (timesp->t1)*1e6);*/		//<== t1 should be displayed

			myovst1234[1] = (timesp->t1)*13.56e6;  /*** changeart= Added save the t1 value to export ***/
			my_ovst1234.t1 = (timesp->t1)*13.56e6;
			/* fprintf(stdout,"t1 = %f/fc  \n", (timesp->t1)*13.56e6);*/
			/* fprintf(stdout,"t2 = %f microsec. \n", (timesp->t2)*1e6);*/		//<== t2 should be displayed

			myovst1234[2] = (timesp->t2)*13.56e6;  /*** changeart= Added save the t2 value to export ***/
			my_ovst1234.t2 = (timesp->t2)*13.56e6;
			/* fprintf(stdout,"t2 = %f/fc  \n", (timesp->t2)*13.56e6);*/
			/*fprintf(stdout,"t3 = %f microsec. \n", (timesp->t3)*1e6);*/		//<== t3 should be displayed

			myovst1234[3] = (timesp->t3)*13.56e6;  /*** changeart= Added save the t3 value to export ***/
			my_ovst1234.t3 = (timesp->t3)*13.56e6;
			/* fprintf(stdout,"t3 = %f/fc  \n", (timesp->t3)*13.56e6);*/
			/* fprintf(stdout,"t4 = %f microsec. \n", (timesp->t4)*1e6);*/		//<== t4 should be displayed

			myovst1234[4] = (timesp->t4)*13.56e6;   /*** changeart= Added save the t3 value to export ***/
			my_ovst1234.t4 = (timesp->t4)*13.56e6;
			/*fprintf(stdout,"t4 = %f/fc  \n", (timesp->t4)*13.56e6);
			fprintf(stdout,"---amplitudes---------\n");
			fprintf(stdout,"Hmax = %f volts \n", Hmax);
			fprintf(stdout,"Max. Amplitude = %f volts \n", (shootreader2->above));*/

		}
				   break;
		case (212) :
		case (424) :
		case (848) :
		{
			/*fprintf(stdout,"hovs = %f \n", (((shootreader2->above-Hmax)/Hmax)));
			fprintf(stdout,"---timings------------\n");
			fprintf(stdout,"t1 = %f microsec. \n", (timesp->t1)*1e6);
			fprintf(stdout,"t1 = %f/fc  \n", (timesp->t1)*13.56e6);
			fprintf(stdout,"t5 = %f microsec. \n", (timesp->t5)*1e6);
			fprintf(stdout,"t5 = %f/fc  \n", (timesp->t5)*13.56e6);
			fprintf(stdout,"t6 = %f microsec. \n", (timesp->t6)*1e6);
			fprintf(stdout,"t6 = %f/fc  \n", (timesp->t6)*13.56e6);
			fprintf(stdout,"---amplitudes---------\n");
			fprintf(stdout,"Hmax = %f volts \n", Hmax);
			fprintf(stdout,"a = %f %% of Hinitial \n", ((timesp->a)/Hmax));*/
		}
				   break;
		}
	}
	break;
	case 'B':
	{
		/* fprintf(stdout,"---RESULTS------------\n");
		fprintf(stdout,"Type B - Bitrate %d\n", rate);
		fprintf(stdout,"---timings------------\n");
		fprintf(stdout,"tf = %f microsec. \n", (timesp->tf)*1e6);
		fprintf(stdout,"tf = %f/fc  \n", (timesp->tf)*13.56e6);
		fprintf(stdout,"tr = %f microsec. \n", (timesp->tr)*1e6);
		fprintf(stdout,"tr = %f/fc  \n", (timesp->tr)*13.56e6);
		fprintf(stdout,"---modulation---------\n");
		fprintf(stdout,"m = %f %% \n", m);
		fprintf(stdout,"---amplitudes---------\n");
		fprintf(stdout,"a = %f volts \n", Hmax);
		fprintf(stdout,"b = %f volts \n", timesp->b);
		fprintf(stdout,"---Overshoots---------\n");
		fprintf(stdout,"hf = %f %% of Hinitial-b\n", (shootreader2->hf_reader)*100);
		fprintf(stdout,"hr = %f %% of Hinitial-b\n", (shootreader2->hr_reader)*100);
		ovsb1=(timesp->b-(shootreader2->above_b))*1000;
		ovsb2=((shootreader2->above)-Hmax)*1000;
		if (ovsb1>0)
		fprintf(stdout,"hf = %f milivolts \n", ovsb1);
		else
		fprintf(stdout,"hf = 0 milivolts \n");
		if (ovsb2>0)
		fprintf(stdout,"hr = %f milivolts \n", ovsb2);
		else
		fprintf(stdout,"hr = 0 milivolts \n");*/
	}
	break;
	}

	return(my_ovst1234);
}

/**************************************************ADD new*********************************************************/
OVST1234 expret(char* csvpath, char type, int rate, double* exportvolts, double* exporttimes, double* ovst1234)
{
	/********************************************ADD NEW*********************************************/

	char voltstr[25];          // intermediate char array to modify the voltage values
	char timestr[25];          // intermediate char array to modify the time values
	double snum = 0;
	double tnum = 0;
	double t = 0;
	int filterlength = 0;
	double Hmax = 0;
	double Hmin = 0;
	double Hmax2 = 0;
	double Hmin2 = 0;
	
	FILE *input_u2 = NULL;
	FILE *poutput = NULL;
	double m = 0.0;
	int length = 0;
	double val = 0;
	int posval = 0;
	int negval = 0;
	double tini = 0;
	double tfin = 0;
	int samples = 0;
	int out_i = 0;
	int length_total = 0;
	int sample_ini = 0;
	int sample_end = 0;
	int flag_cut = 0;
	int samplesp = 0;
	int fi = 0;                 // Filter generic index
	double b1 = 0;            // Filter parameters
	double b2 = 0;
	double b3 = 0;
	double b4 = 0;
	double b5 = 0;
	double a1 = 0;
	double a2 = 0;
	double a3 = 0;
	double a4 = 0;
	double a5 = 0;
	double freq1 = 0;
	double freq2 = 0;
	double as[5] = { 0 };
	double bs[5] = { 0 };
	double t0 = 0;
	double tlast = 0;
	int lineskip = 0;
	char* infomsg;  /**Added saved the printf messages **/
	
	double vv_ovs = 0.0;
	double vv_t1 = 0.0;
	double vv_t2 = 0.0;
	double vv_t3 = 0.0;
	double vv_t4 = 0.0;
	OVST1234 mmyovst1234;
	mmyovst1234.ovs = 0;
	mmyovst1234.t1 = 0;
	mmyovst1234.t2 = 0;
	mmyovst1234.t3 = 0;
	mmyovst1234.t4 = 0;
	*exportvolts = 0;
	*exporttimes = 0;
	
	ovst1234 = NULL;
	/********************************************ADD NEW***********************************************/
	int j = 0;
	float timevalue = 0;
	float voltvalue = 0;
	
	char volt[25];
	char time[25];
	double volt1[20000];
	double time1[20000];
	double myovst1234[5] = { 0.0, 0.0, 0.0, 0.0, 0.0 };
	
	infomsg = "";
	FILE *pointfile = NULL;
	if (type == 'B' || rate == 212 || rate == 424 || rate == 848)
		infomsg = "this case is not implemented";

	

	if (type == 'A' && rate == 106)
	{ 
		
		pointfile = fopen(csvpath, "r");
		for (int j = 0; j < 10; j++)
			skip_line(pointfile);


		while (!feof(pointfile))                  // read the input file
		{
			read_line(pointfile, volt, time);
			exportvolts[j] = atof(volt);
			exporttimes[j] = atof(time);
			j++;
		}
		fclose(pointfile);
	}
	/***************************************ADD NEW***********************************************/
	double *voutput = malloc(sizeof(double)*MAX_SAMPLES);
	double *toutput = malloc(sizeof(double)*MAX_SAMPLES);
	double *envelope = malloc(sizeof(double)*MAX_SAMPLES);
	double *vfilter = malloc(sizeof(double)*MAX_SAMPLES);
	double *tfilter = malloc(sizeof(double)*MAX_SAMPLES);
	TIMES *timesp = (TIMES *)malloc(sizeof(TIMES));
	TIMES *timesp2 = (TIMES *)malloc(sizeof(TIMES));
	SHOOTREADER *shootreader2 = (SHOOTREADER *)malloc(sizeof(SHOOTREADER));

	if (voutput != NULL && toutput != NULL && envelope != NULL && vfilter != NULL && tfilter != NULL && timesp != NULL && timesp2 != NULL && shootreader2 != NULL)
	{
		memset(voutput, 0, MAX_SAMPLES);
		memset(toutput, 0, MAX_SAMPLES);
		memset(envelope, 0, MAX_SAMPLES);
		memset(vfilter, 0, MAX_SAMPLES);
		memset(tfilter, 0, MAX_SAMPLES);

		if ((type != 'A' && type != 'B') || (rate != 106 && rate != 212 && rate != 424 && rate != 848))
		{
			//fprintf(stdout, "Wrong Type (A or B) or Bitrate (106, 212, 424, 848)");
			infomsg = "Wrong Type (A or B) or Bitrate (106, 212, 424, 848)";
		}
		else
		{
			/*pointfile=fopen(argv[3],"r");*/
			pointfile = fopen(csvpath, "r");
			input_u2 = fopen("intermediate.txt", "w");            // modified-intermediate amplitude vector

			if (pointfile != NULL && input_u2 != NULL)
			{
				//1. LOAD DATA + CHECKING DATA (WITHOUT FILTER)
				for (lineskip = 0; lineskip<10; lineskip++)      // Skips the first 10 lines which are the header of csv files
				{
					skip_line(pointfile);
				}
				read_line(pointfile, voltstr, timestr);
				t0 = atof(timestr);
				while (!feof(pointfile))                  // We are reading the lines of the voltage input file
				{
					if (voltstr[0] != '\0')
					{
						snum = atof(voltstr);
						tnum = atof(timestr);
						if (snum<0)
							negval++;
						else
							posval++;
						vfilter[samplesp] = snum;
						tfilter[samplesp] = tnum;
						samplesp++;
						read_line(pointfile, voltstr, timestr);
					}
					tlast = tfilter[samplesp - 1];
				}
				samplesp = samplesp + 3;

				samplesp = datacheck(posval, negval, samplesp, tlast, pointfile);

				tlast = tfilter[samplesp];

				//2. DATA FILTER 10 MHz BANDPASS
				//butter_order=2;  // Real 4th Order (if !=2, the function won큧 work!!!!)
				freq1 = 8.56e6 / (1 / (2 * ((tlast - t0) / (samplesp - 1))));
				freq2 = 18.56e6 / (1 / (2 * ((tlast - t0) / (samplesp - 1))));
				butterworth_coeffs(freq1, freq2, as, bs);

				b1 = bs[0];
				b2 = bs[1];
				b3 = bs[2];
				b4 = bs[3];
				b5 = bs[4];

				a1 = as[0];
				a2 = as[1];
				a3 = as[2];
				a4 = as[3];
				a5 = as[4];

				for (fi = 0; fi<samplesp; fi++)
				{
					if (fi<7 || fi>samplesp - 7)
						voutput[fi] = 0;
					else
						voutput[fi] = (b1*vfilter[fi] + b2*vfilter[fi - 1] + b3*vfilter[fi - 2] + b4*vfilter[fi - 3] + b5*vfilter[fi - 4] - a2*voutput[fi - 1] - a3*voutput[fi - 2] - a4*voutput[fi - 3] - a5*voutput[fi - 4]) / a1;
				}

				rewind(pointfile);
				lineskip = 0;
				for (lineskip = 0; lineskip<10; lineskip++)      // Skips the first 10 lines which are the header of csv files
				{
					skip_line(pointfile);
				}
				for (fi = 0; fi<(samplesp - 7); fi++)       /* We are reading the lines of the voltage input file */
				{
					val = voutput[fi];
					read_line(pointfile, voltstr, timestr);
					fprintf(input_u2, "%s,%f\n", timestr, val);
					length++;
				}

				//3. HILBERT TRANSFORM AND THE COMPLEX ENVELOPE
				rewind(input_u2);
				hilbert("intermediate.txt");         // performs hilbert transform

				poutput = fopen("output.txt", "r");       // hilbert transform output vector

				read_line(poutput, voltstr, timestr);
				tini = atof(timestr);
				rewind(poutput);

				if (poutput != NULL)
				{
					while (!feof(poutput))    // We are reading the lines of the voltage input file */
					{
						read_line(poutput, voltstr, timestr);
						if (timestr[0] != '\0')
						{
							snum = atof(voltstr);
							voutput[samples] = snum;
							t = atof(timestr);
							toutput[samples] = t;
							samples++;//==>US  // Same variable as the one in Hmaxfinder
							tfin = t;
						}
					}
				}
				else
					infomsg = "Error in Hilbert transform\n";

				//4. USING A SMOOTHING FILTER (MOV. AVG) TO REDUCE THE NOISE
				filterlength = 1;
				length_total = envfilt(voutput, toutput, filterlength, tini, tfin, samples, envelope);

				//5. 100% OF H_INITIAL
				Hmaxfinder(envelope, &Hmax, &Hmin, length_total);

				//6. COMPUTING THE ISO BASED TIMES
				tfinder(type, envelope, toutput, tini, Hmax, Hmin, rate, length_total, timesp);

				//7. CHECKING FOR ISO DEFINED MONOTONY
				monocheck(envelope, toutput, Hmax, timesp, rate, type);

				out_i = 0;
				while (out_i<MAX_SAMPLES)      // Finds how many zeros are at the beginning of vector envelope
				{
					if (envelope[out_i] == 0 && flag_cut == 0)
					{
						sample_ini = out_i;
						tini = toutput[sample_ini + 1];
					}

					if (envelope[out_i] != 0)
					{
						flag_cut = 1;
						sample_end = out_i;
						tfin = toutput[sample_end];
					}
					out_i++;
				}

				samples = sample_end - sample_ini - 1;      //==>US

				for (out_i = 0; out_i<samples; out_i++)
				{
					voutput[out_i] = envelope[out_i + sample_ini + 1];
					toutput[out_i] = toutput[out_i + sample_ini + 1];
				}
				for (out_i = samples + 1; out_i<MAX_SAMPLES; out_i++)
				{
					voutput[out_i] = 0.0;
					toutput[out_i] = 0.0;
				}

				tini = toutput[0];
				tfin = toutput[samples];

				//8. OVERSHOOT OF THE READER
				//  fprintf (stdout,"\n");    // 2nd set of functions, "New Line" printed for debug purposes
				filterlength = 3;
				length_total = envfilt(voutput, toutput, filterlength, tini, tfin, samples, envelope);      // 2nd Filtering to find the alternate envelope
				Hmaxfinder(envelope, &Hmax2, &Hmin2, length_total);
				tfinder(type, envelope, toutput, tini, Hmax2, Hmin2, rate, length_total, timesp2);
				monocheck(envelope, toutput, Hmax2, timesp2, rate, type);               // The parameters of the alternate envelope are calculated
				overshoot(timesp2, Hmax2, envelope, toutput, rate, type, samples, shootreader2);   // This time the over- and undershoots are found

				//9. ISO RESTRICTION FOR THE SHOOT
				// (NOT IMPLEMENTED)

				//10. MODULATION
				m = modulation(type, Hmax, timesp->b);

				//11. DISPLAY
				mmyovst1234 = display(type, rate, shootreader2, timesp, Hmax, m, &myovst1234);
				mmyovst1234.msg = infomsg;
				ovst1234 = &myovst1234;
				
			}
			else if (pointfile == NULL || input_u2 != NULL)
				fprintf(stdout, "file(s) could not be opened \n");

			fclose(pointfile);
			fclose(input_u2);
		}
	}

	else
		infomsg = "Memory could not be allocated\n"; //fprintf(stdout, "Memory could not be allocated");

	free(voutput);
	free(toutput);
	free(envelope);
	free(vfilter);
	free(tfilter);
	free(timesp);
	free(timesp2);
	free(shootreader2);

	/***************************************ADD NEW***********************************************/
	
	return(mmyovst1234);
}
