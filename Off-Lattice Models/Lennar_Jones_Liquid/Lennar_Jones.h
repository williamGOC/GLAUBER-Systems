#ifndef LENNAR_JONES_H
#define LENNAR_JONES_H

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <assert.h>

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

#define L 27.386127875258
#define NC 15
#define RCUT 2.5
#define N 225

#define DMAX 0.5
#define BETA 0.8

typedef double vectorR[2];
typedef vectorR *vectorRN;

typedef struct particles{
   vectorRN X;
   int step;
   int Naccep;
   long seed;
   double E;
} world;

world *makeSystem();
void freeWorld(world *pW);
void printerSystem(world *pW, FILE *gnuplot);


void pbc(world *pW);
double square(double x);
double minImage(double xi, double xj);
double distij(vectorR ri, vectorR rj);
double energy(world *pW);

void stepMC(world *pW);
float rand2(long *idum);
#endif // LENNAR_JONES_H