#include "ising_cpu.h"


IsingModel makeSystem(double _J, double _T, double _H){

    // Memory allocation for the central structure of the program.
	IsingModel pS = (IsingModel)malloc(sizeof(struct system));		
	assert(pS);

	pS -> N = (int)pow(L, DIM);

	// Creation of dynamic array member of the central structure.
	pS -> spins = (int *)malloc(pS -> N * sizeof(int));
	assert(pS -> spins);

	pS -> seed = -time(NULL);
	pS -> step = 0;

	pS -> J = _J;
	pS -> T = _T;
	pS -> H = _H;

	
	// Calculate lattice volume elements:
	unsigned volume[DIM+1];
	unsigned i; for (i = 0; i <= DIM; ++i){
		if (i == 0)
			volume[i] = 1;
		else
			volume[i] = volume[i-1] * L;
	}


	// Neighborhood table:
	pS -> neighbors = (unsigned *)malloc(2 * DIM * pS -> N *sizeof(unsigned));
	assert(pS -> neighbors);

	
	
	for (i = 0; i < pS -> N; ++i){
		unsigned short k = 0;

		unsigned short dim; for (dim = 0; dim < DIM; ++dim){					// dimension loop
			short dir; for (dir = -1; dir <= 1; dir += 2){ 			            // two directions in each dimension

				// neighbor's id in spin list:
				int npos = i + dir * volume[dim];

				// periodic bonduary conditions:
				int test = (i % volume[dim+1]) + dir * volume[dim];

				if (test < 0)
					npos += volume[dim+1];
				else if (test >= volume[dim+1])
					npos -= volume[dim+1];

				pS -> neighbors[2*DIM*i + k] = npos;
				k++;
			}
		}

		pS -> spins[i] =  (int)(2 * floor(2 * rand2(&pS-> seed)) - 1);   // create a random system
		//pS -> spins[i] =  UP;   // create a UP-system	
	}

	return pS;
}



void stepMC(IsingModel pS){


	unsigned i, j; for (i = 0; i < pS -> N; ++i){
		
		unsigned index = (unsigned)(pS -> N * rand2(&pS -> seed));
		int spinstate = pS -> spins[index];

		double deltaE = 0.0;

		for (j = 0; j < 2 * DIM; ++j)
			deltaE += pS -> J * pS -> spins[pS -> neighbors[2*DIM*index + j]];

		deltaE = (pS -> step > 0.2*N_SWEEPS)?(2.0 * spinstate * (deltaE + pS -> H * cos((2 * M_PI * (pS -> step - 0.2*N_SWEEPS)) / pS -> T))):(2.0 * spinstate * deltaE);
		

		if (exp(-1.0 * deltaE) >= rand2(&pS -> seed))
			pS -> spins[index] = - spinstate;

	}

	++pS -> step;
}

double getm(IsingModel pS){

	int M = 0;
	int i; for (i = 0; i < pS -> N; ++i)
		M += pS -> spins[i];

	return (1.0 * M) / (1.0 * pS -> N);
}


void destroySystem(IsingModel pS){

	free(pS -> spins);
	free(pS -> neighbors);
	free(pS);
}


unsigned getId(int x, int y){
	return x * L + y;
}



void printerSystem(IsingModel pS, FILE * gnuplotPIPE){

	fprintf(gnuplotPIPE, "set title '{/=20 Ising Model}'\n");
	fprintf(gnuplotPIPE, "unset xtics; unset ytics\n");
	fprintf(gnuplotPIPE, "unset border\n");
	fprintf(gnuplotPIPE, "set xrange [0:%d]; set yrange [0:%d]\n", L, L);
	fprintf(gnuplotPIPE, "plot '-' u 1:2:3 w p pt 5 ps 3 lc variable title ''\n");

	int i; for(i = 0; i < pS -> N; i++)
		fprintf(gnuplotPIPE, "%d\t%d\t%d\n", i%L, i/L, pS -> spins[i]);

	fprintf(gnuplotPIPE, "e\n");
	fflush(gnuplotPIPE);
}

void addSquare(int x, int y, int w, int sign, IsingModel pS){
		
	int i, j;
	int ww=(w>L)?(L):(w);
		
	for(j=y-ww/2;j<y+ww/2;j++){
		for(i=x-ww/2;i<x+ww/2;i++){
			pS -> spins[(i+L)%L+((j+L)%L)*L]=sign;
		}
	}
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
