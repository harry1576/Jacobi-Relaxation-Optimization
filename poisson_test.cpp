
/// \brief Test program for solving Poisson's equation using Jacobi relaxation
/// \author M. P. Hayes UCECE

#include <stdlib.h>
#include <stdio.h>
#include "poisson.hpp"
#include <pthread.h>
#include <stdint.h>
#include <inttypes.h>

int main (int argc, char *argv[])
{
    double *source;
    double *potential;
    unsigned int N;
    unsigned int numiters;
    unsigned int numcores;    
    unsigned int xsize;
    unsigned int ysize;
    unsigned int zsize;    
    double delta = 0.1;
    uint64_t nanoseconds;
    struct timespec start, end;
	size_t count = 1000000000;


    if (argc < 3)
    {
        fprintf (stderr, "Usage: %s size numiters\n", argv[0]);
        return 1;
    }

    N = atoi(argv[1]);
    xsize = N;
    ysize = N;
    zsize = N;

    numiters = atoi(argv[2]);

    if (argc > 3)
        numcores = atoi(argv[3]);
    else
        numcores = 0;

    source = (double *)calloc(xsize * ysize * zsize, sizeof(*source));
    potential = (double *)calloc(xsize * ysize * zsize, sizeof(*potential));

    source[((zsize / 2 * ysize) + ysize / 2) * xsize + xsize / 2] = 1.0;    
    
    clock_gettime(CLOCK_MONOTONIC, &start);
    
    poisson_dirichlet(source, potential, 0, xsize, ysize, zsize, delta,numiters, numcores);
    
    
   //free(source);
   //free(potential);

    
    
    clock_gettime(CLOCK_MONOTONIC, &end);
    printf("Took %" PRIu64 " ms for %zu iterations at size %d and cores %zu \n",nanoseconds / 1000000, numiters,numcores);

	nanoseconds = (end.tv_sec - start.tv_sec) * 1000000000ULL + (end.tv_nsec - start.tv_nsec);
	
	
	

    
	//for (int i = 0; i < xsize*ysize*zsize; i++) {             
    //for (int z = 0; z < zsize; ++z) {
		//for (int y = 0; y < ysize; ++y) {
			//for (int x = 0; x < xsize; ++x) {
				//printf("%f ", potential[((z * ysize) + y) * xsize + x]);
			//}
			//printf("\n");
		//}
		//printf("\n\n");
	//}
	
	//free(source);
	//free(potential);


    

    return 0;
}
