#include "glauber_system.c"


int main(int argc, char **argv){

	PottsModel pP;

	pP = makeSystem(2, 0.48);

	
	srand(time(NULL));


	printf("%lf\n", pP -> Q * rand2(&pP-> seed));
	/*int i; for (i = 0; i < pP -> N; ++i){
			printf("%d\n", pP -> spins[i]);

	}*/
	
	
	destroySystem(pP);

	return 0;
}