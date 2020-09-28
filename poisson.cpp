#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
#include <thread>
#include <vector>
#include <time.h>

void poisson_dirichlet_test(double * __restrict__ source,
							double * __restrict__ potential,
							double Vbound,
							unsigned int xsize, unsigned int ysize, unsigned int zsize, double delta,
							unsigned int numiters, unsigned int numcores)
{
	// source[i, j, k] is accessed with source[((k * ysize) + j) * xsize + i]
    // potential[i, j, k] is accessed with potential[((k * ysize) + j) * xsize + i]    
    size_t size = (size_t)ysize * zsize * xsize * sizeof(double);
	double *input = (double *)malloc(size);
	if (!input) {
		fprintf(stderr, "malloc failure\n");
		return;
	}
	memcpy(input, source, size);
	for (unsigned int iter = 0; iter < numiters; iter++) {
		for (unsigned int x = 0; x < xsize; x++) {
			for (unsigned int z = 0; z < zsize; z++) {
				for (unsigned int y = 0; y < ysize; y++) {
					double res = 0;

					if (x < xsize - 1)
						res += input[((z * ysize) + y) * xsize + (x + 1)];
					if (x > 0)
						res += input[((z * ysize) + y) * xsize + (x - 1)];

					if (y < ysize - 1)
						res += input[((z * ysize) + (y + 1)) * xsize + x];
					if (y > 0)
						res += input[((z * ysize) + (y - 1)) * xsize + x];

					if (z < zsize - 1)
						res += input[(((z + 1) * ysize) + y) * xsize + x];
					if (z > 0)
						res += input[(((z - 1) * ysize) + y) * xsize + x];

					res -= delta * delta * source[((z * ysize) + y) * xsize + x];

					res /= 6;

					potential[((z * ysize) + y) * xsize + x] = res;
				}
			}
		}
		memcpy(input, potential, size);
	}
	free(input);
}


int calculate_index(int i, int j ,int k, int xsize, int ysize)
{
	return (((k * ysize) + j) * xsize + i);
}


void poisson_thread_function(int startIndex, int endIndex)
{
	
}



/// Solve Poisson's equation for a rectangular box with Dirichlet
/// boundary conditions on each face.
/// \param source is a pointer to a flattened 3-D array for the source function
/// \param potential is a pointer to a flattened 3-D array for the calculated potential
/// \param Vbound is the potential on the boundary
/// \param xsize is the number of elements in the x-direction
/// \param ysize is the number of elements in the y-direction
/// \param zsize is the number of elements in the z-direction
/// \param delta is the voxel spacing in all directions
/// \param numiters is the number of iterations to perform
/// \param numcores is the number of CPU cores to use.  If 0, an optimal number is chosen
void poisson_dirichlet (double * __restrict__ source,
                        double * __restrict__ potential,
                        double Vbound,
                        unsigned int xsize, unsigned int ysize, unsigned int zsize, double delta,
                        unsigned int numiters, unsigned int numcores)
{
    // Need equal size threads for now.
	//unsigned int size = xsize*ysize*zsize;
	//if (size % numcores) {s
	//	fprintf(stderr, "size %zu must be divisible by %d\n", size, numcores);
	//	return 0;
	//}
	
	// Array of threads
	//thread threads[numcores];
	
	// Split up data and spawn the threads
	//for (int i = 0; i < numcores; i++) {
		//threads[i] = thread(poisson_dirichlet_func, *
	//poisson_dirichlet_func(* __restrict__ source,  * __restrict__ potential, 0, xsize, ysize, zsize, delta, numiters, numcores);
	
	
	
	std::vector<std::thread> threads;
	int startIndex; //Start index for each core...
	int endIndex; // End index for each core...
	
	
	// Generates threads equal to the number of cores...
	for (unsigned int i = 0; i < numcores; ++i) {
		startIndex = (i * xsize / numcores) + 1;
		endIndex = (i + 1) * xsize / numcores + 1;
		threads.push_back(std::thread(poisson_thread_function,startIndex,endIndex));
	}
	
	
	
	
	
	
	// Join threads...
	for (auto& thread : threads) {
			thread.join();
	}
	

	
	
	
	
}


