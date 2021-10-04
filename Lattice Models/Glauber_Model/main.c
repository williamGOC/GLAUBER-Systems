#include "glauber_system.c"


int main(int argc, char **argv){

	GlauberModel pG;

	pG = makeSystem(0.48, 1000, 0.4);

	
	srand(time(NULL));


	
	int i; for (i = 0; i < N_SWEEPS; ++i){
		stepMC(pG);
		
		if (i > 0.2 * N_SWEEPS)
			printf("%d\t%lf\n", i - (int)(0.2*N_SWEEPS), getm(pG));

	}
	
	
	destroySystem(pG);

	return 0;
}