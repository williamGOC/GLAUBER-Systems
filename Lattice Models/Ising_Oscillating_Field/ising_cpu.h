#ifndef _ising_cpu_h_
#define _ising_cpu_h_

#include <time.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#define UP 1
#define DOWN -1

#define DIM 2
#define L 500

#define N_SWEEPS 10000

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

typedef struct system{
	
	int *spins;
	unsigned *neighbors;

	int N;
	long seed;
	int step;

	double J, T, H;

} *IsingModel;


// System Structuring Functions:
IsingModel makeSystem(double _J, double _T, double _H);
void destroySystem(IsingModel pS);
void printerSystem(IsingModel pS, FILE * gnuplotPIPE);


// Dynamical Evolution Functions:
void stepMC(IsingModel pS);
float rand2(long *idum);
unsigned getId(int x, int y);
double getm(IsingModel pS);
void addSquare(int x, int y, int w, int sign, IsingModel pS);
#endif  // _ising_cpu_h_