#include <time.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>


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

	int N, Q;
	long seed;

	double J;
} *PottsModel;


// System Structuring Functions:
PottsModel makeSystem(double _Q, double _J);
void destroySystem(PottsModel pP);

// Dynamical Evolution Functions:
void stepMC(PottsModel pP);
float rand2(long *idum);
unsigned getId(int x, int y);



PottsModel makeSystem(double _Q, double _J){


	// Memory allocation for the central structure of the program.
	PottsModel pP = (PottsModel)malloc(sizeof(struct system));		
	assert(pP);

	pP -> N = (int)pow(L, DIM);

	// Creation of dynamic array member of the central structure.
	pP -> spins = (int *)malloc(pP -> N * sizeof(int));
	assert(pP -> spins);

	pP -> seed = -time(NULL);

	pP -> J = _J;
	pP -> Q = _Q;

	
	// Calculate lattice volume elements:
	unsigned volume[DIM+1];
	unsigned i; for (i = 0; i <= DIM; ++i){
		if (i == 0)
			volume[i] = 1;
		else
			volume[i] = volume[i-1] * L;
	}


	// Neighborhood table:
	pP -> neighbors = (unsigned *)malloc(2 * DIM * pP -> N *sizeof(unsigned));
	assert(pP -> neighbors);

	
	
	for (i = 0; i < pP -> N; ++i){
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

				pP -> neighbors[2*DIM*i + k] = npos;
				k++;
			}
		}
			
		pP -> spins[i] = floor(pP -> Q * rand2(&pP-> seed));   // create a random system	
	}

	return pP;
}



void stepMC(PottsModel pP){


	unsigned i, j; for (i = 0; i < pP -> N; ++i){
		
		unsigned index = (unsigned)(pP -> N * rand2(&pP -> seed));
		int spinstate  = pP -> spins[index];
		int newstate   = floor(pP -> Q * rand2(&pP-> seed));

		double Ebefore = 0.0;
		double Eafter  = 0.0;

		for (j = 0; j < 2 * DIM; ++j){
			if (spinstate == pP -> spins[pP -> neighbors[2 * DIM * index + j]])
				Ebefore += pP -> J;
			if (newstate == pP -> spins[pP -> neighbors[2 * DIM * index + j]])
				Eafter += pP -> J;
		}
			
		if (exp(-1.0 * (Eafter - Ebefore)) >= rand2(&pP -> seed))
			pP -> spins[index] = newstate;

	}
}



void destroySystem(PottsModel pP){

	free(pP -> spins);
	free(pP -> neighbors);
	free(pP);
}


unsigned getId(int x, int y){
	return x * L + y;
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