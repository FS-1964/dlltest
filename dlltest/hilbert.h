/***************************************************************************/
/*hilbert.h                                                                */
/*This code contains the necessary functions for extracting envelope       */
/***************************************************************************/

#ifndef HILBERT_H_
#define HILBERT_H_

/*This function reads the sampled data recorded in the file*/
int ReadData(void);

/*This function function performs the fourier transform*/
void Fft(void);

/*This function performs the necessary phase shift*/
void PhaseShifting(void);

/*This function performs the inverse fourier transfor*/
void Ifft(void);

/*Envelope reconstruction is done by this function*/
int EnvelopeReconstruction(void);

/*Hilbert main function*/
void hilbert(char *fnamep);

#endif /* HILBERT_H_ */

