#include "Lennar_Jones.h"

world *makeSystem(){

    world *pW = (world *)malloc(sizeof(world));
    assert(pW);

    pW -> X = (vectorRN)malloc(N * sizeof(vectorR));
    assert(pW -> X);

    pW -> step = 0;
    pW -> seed = -time(NULL);
    pW -> Naccep = 0;

    int i; for (i = 0; i < N; i++){
        pW -> X[i][0] = (i % NC) * L / NC;
        pW -> X[i][1] = (i / NC) * L / NC;
    }

    pW -> E = 0.0;
    double dij, tem;

    int j; for (i = 0; i < N; i++){
        for (j = i + 1; j < N; j++){
            dij = distij(pW -> X[i], pW -> X[j]);
            tem = pow(dij, -6);
            pW -> E += 4 * tem * (tem - 1);
        }
    }
    return pW;
}

void freeWorld(world *pW){

    free(pW -> X);
    free(pW);
}


void pbc(world *pW){

    double size = 1.0 / L;
    int i; for (i = 0; i < N; i++){
        pW -> X[i][0] -= floor(pW -> X[i][0] * size) * L;
        pW -> X[i][1] -= floor(pW -> X[i][1] * size) * L;
    }
}


double square(double x){
    return x * x;
}

double mimImage(double xi, double xj){
    double xij = xi - xj;
    return xij - L * round(xij / L);
}

double distij(vectorR ri, vectorR rj){
    return sqrt(square(mimImage(ri[0], rj[0])) + square(mimImage(ri[1], rj[1])));
}

double energy(world *pW){

    double U = 0.0;
    double dij, tem;

    int i, j; for (i = 0; i < N; i++){
        for (j = i + 1; j < N; j++){
            dij = distij(pW -> X[i], pW -> X[j]);
            tem = pow(dij, -6);
            U += 4 * tem * (tem - 1);
        }
    }
    
    return U;
}


void stepMC(world *pW){
        
    unsigned index = (unsigned)(N * rand2(&pW -> seed));
    double Vold = pW -> E;

    float dir0 = 2.0 * rand2(&pW -> seed) - 1.0;
    float dir1 = 2.0 * rand2(&pW -> seed) - 1.0;

    pW -> X[index][0] += DMAX * dir0;
    pW -> X[index][1] += DMAX * dir1;


    double Vnew = energy(pW);
    double DELTAV = Vnew - Vold;
    double DELTAVB = BETA * DELTAV;

    if (DELTAVB < 75.0){
        if (DELTAVB <= 0.0){
            pW -> Naccep ++;
            pW -> E += DELTAV;
        }
        else if (exp(-1.0 * DELTAVB) > rand2(&pW -> seed)){
            pW -> Naccep ++;
            pW -> E += DELTAV;
        }
        else {
            pW -> X[index][0] -= DMAX * dir0;
            pW -> X[index][1] -= DMAX * dir1;
        }
    }
    else{
        pW -> X[index][0] -= DMAX * dir0;
        pW -> X[index][1] -= DMAX * dir1;
    } 

    pW -> step ++;
}


void printerSystem(world *pW, FILE *gnuplot){
   
    fprintf(gnuplot, "set title '{/=20 Gas Model, pass %d}'\n", pW -> step);
    fprintf(gnuplot, "set xlabel  '{/=15 X}'\n");  
    fprintf(gnuplot, "set ylabel  '{/=15 Y}'\n");

    fprintf(gnuplot, "set xrange [0:%lf]\n", L);
    fprintf(gnuplot, "set yrange [0:%lf]\n", L);
    fprintf(gnuplot, "plot '-' with p pointtype 6 t ''\n");
  
    int i; for(i = 0; i < N; i++)
        fprintf(gnuplot, "%lf\t%lf\n", pW -> X[i][0], pW -> X[i][1]);
    
    fprintf(gnuplot, "e\n");
    fflush(gnuplot);
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