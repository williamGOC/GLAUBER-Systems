#include "Lennar_Jones.h"

int main(int argc, char **argv) {

    world *pW = makeSystem();
    srand(time(NULL));

    FILE *pipe = popen("gnuplot", "w");
    assert(pipe);

    int i; for (i = 0; i < 100000; i++){
        printerSystem(pW, pipe);
        stepMC(pW);
        pbc(pW);
    }

    freeWorld(pW);
    fclose(pipe);

    return 0;
}
