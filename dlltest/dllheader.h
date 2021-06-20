#ifdef DLLHEADER_EXPORTS
#define DLLHEADER_API __declspec(dllexport) 
#else
#define DLLHEADER_API __declspec(dllimport) 
#endif

#include "ovst1234.h"



DLLHEADER_API OVST1234 expret(char* csvpath, char typ, int rate, double* exportvolts, double* exporttimes, double* calculatedvalues);