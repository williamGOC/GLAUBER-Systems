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

} *GlauberModel;


// System Structuring Functions:
GlauberModel makeSystem(double _J, double _omeg, double _H);
void destroySystem(GlauberModel pG);
void printerSystem(GlauberModel pG, FILE * gnuplotPIPE);


// Dynamical Evolution Functions:
void stepMC(GlauberModel pG);
float rand2(long *idum);
unsigned getId(int x, int y);
double getm(GlauberModel pG);



GlauberModel makeSystem(double _J, double _T, double _H){


	// Memory allocation for the central structure of the program.
	GlauberModel pG = (GlauberModel)malloc(sizeof(struct system));		
	assert(pG);

	pG -> N = (int)pow(L, DIM);

	// Creation of dynamic array member of the central structure.
	pG -> spins = (int *)malloc(pG -> N * sizeof(int));
	assert(pG -> spins);

	pG -> seed = -time(NULL);
	pG -> step = 0;

	pG -> J = _J;
	pG -> T = _T;
	pG -> H = _H;

	
	// Calculate lattice volume elements:
	unsigned volume[DIM+1];
	unsigned i; for (i = 0; i <= DIM; ++i){
		if (i == 0)
			volume[i] = 1;
		else
			volume[i] = volume[i-1] * L;
	}


	// Neighborhood table:
	pG -> neighbors = (unsigned *)malloc(2 * DIM * pG -> N *sizeof(unsigned));
	assert(pG -> neighbors);

	
	
	for (i = 0; i < pG -> N; ++i){
		unsigned short k = 0;

		unsigned short dim; for (dim = 0; dim < DIM; ++dim){					// dimension loop
			short dir; for (dir = -1; dir <= 1; dir += 2){ 			// two directions in each dimension

				// neighbor's id in spin list:
				int npos = i + dir * volume[dim];

				// periodic bonduary conditions:
				int test = (i % volume[dim+1]) + dir * volume[dim];

				if (test < 0)
					npos += volume[dim+1];
				else if (test >= volume[dim+1])
					npos -= volume[dim+1];

				pG -> neighbors[2*DIM*i + k] = npos;
				k++;
			}
		}

		
			
		pG -> spins[i] =  (int)(2 * floor(2 * rand2(&pG-> seed)) - 1);   // create a random system
		//pG -> spins[i] =  UP;   // create a UP-system	
	}

	return pG;
}



void stepMC(GlauberModel pG){


	unsigned i, j; for (i = 0; i < pG -> N; ++i){
		
		unsigned index = (unsigned)(pG -> N * rand2(&pG -> seed));
		int spinstate = pG -> spins[index];

		double deltaE = 0.0;

		for (j = 0; j < 2 * DIM; ++j)
			deltaE += pG -> J * pG -> spins[pG -> neighbors[2*DIM*index + j]];

		deltaE = (pG -> step > 0.2*N_SWEEPS)?(2.0 * spinstate * (deltaE + pG -> H * cos((2 * M_PI * (pG -> step - 0.2*N_SWEEPS)) / pG -> T))):(2.0 * spinstate * deltaE);
		

		if (exp(-1.0 * deltaE) >= rand2(&pG -> seed))
			pG -> spins[index] = - spinstate;

	}

	++pG -> step;
}

double getm(GlauberModel pG){

	int M = 0;
	int i; for (i = 0; i < pG -> N; ++i)
		M += pG -> spins[i];

	return (1.0 * M) / (1.0 * pG -> N);
}


void destroySystem(GlauberModel pG){

	free(pG -> spins);
	free(pG -> neighbors);
	free(pG);
}


unsigned getId(int x, int y){
	return x * L + y;
}



void printerSystem(GlauberModel pG, FILE * gnuplotPIPE){

	fprintf(gnuplotPIPE, "set title '{/=20 Glauber Model}'\n");
	fprintf(gnuplotPIPE, "unset xtics; unset ytics\n");
	fprintf(gnuplotPIPE, "unset border\n");
	fprintf(gnuplotPIPE, "set xrange [0:%d]; set yrange [0:%d]\n", L, L);
	fprintf(gnuplotPIPE, "plot '-' u 1:2:3 w p pt 5 ps 3 lc variable title ''\n");

	int i; for(i = 0; i < pG -> N; i++)
		fprintf(gnuplotPIPE, "%d\t%d\t%d\n", i%L, i/L, pG -> spins[i]);

	fprintf(gnuplotPIPE, "e\n");
	fflush(gnuplotPIPE);
}



float rand2(long *idum){
    
  int j;
  long k;
  
  static long idum2=123456789;
  static long iy=0;
  static long iv[NTAB];
  float temp;

  if (*idum <= 0){                        //Initialize.
  if (-(*idum) < 1) *idum=1;              //Be sure to prevent idum = 0.
  else *idum = -(*idum);
  idum2=(*idum);
  for (j=NTAB+7;j>=0;j--){                //Load the shuffle table (after 8 warm-ups).
  k=(*idum)/IQ1;
  *idum=IA1*(*idum-k*IQ1)-k*IR1;
  if (*idum < 0) *idum += IM1;
  if (j < NTAB) iv[j] = *idum;
  }
  iy=iv[0];
  }
  k=(*idum)/IQ1;                           //Start here when not initializing.
  *idum=IA1*(*idum-k*IQ1)-k*IR1;           //Compute idum=(IA1*idum) % IM1 without
  if (*idum < 0) *idum += IM1;             //overflows by Schrage’s method.
  k=idum2/IQ2;
  idum2=IA2*(idum2-k*IQ2)-k*IR2;           //Compute idum2=(IA2*idum) % IM2 likewise.
  if (idum2 < 0) idum2 += IM2;
  j=iy/NDIV;                               //Will be in the range 0..NTAB-1.
  iy=iv[j]-idum2;                          //Here idum is shuffled, idum and idum2 are
  iv[j] = *idum;                           //combined to generate output.
  if (iy < 1) iy += IMM1;
  if ((temp=AM*iy) > RNMX) return RNMX;    //Because users don’t expect endpoint values.
  else return temp;
}