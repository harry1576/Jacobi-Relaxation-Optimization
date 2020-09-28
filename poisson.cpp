#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
#include <thread>
#include <vector>
#include <time.h>

void poisson_dirichlet_t(double * __restrict__ source,
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


int get_index(int i, int j ,int k, int xsize, int ysize)
{
	return (((k * ysize) + j) * xsize + i);
}


void poisson_thread_function(int startIndex, int endIndex,double * __restrict__ source, double * __restrict__ potential, unsigned int xsize, unsigned int ysize, unsigned int zsize, double delta, double*  __restrict__ result)
{
	double v1,v2,v3,v4,v5,v6,deltaSquared;
	int k,i,j,currentIndex;
	
	deltaSquared = delta * delta;
	
	for(k=startIndex; k < endIndex; k++)
	{
		for(j=1; j < ysize - 1; j++)
		{
			for(i = 1; i < xsize -1; i++)
			{
				v1 = potential[get_index(i+1,j,k,xsize,ysize)];
				v2 = potential[get_index(i-1,j,k,xsize,ysize)];
				v3 = potential[get_index(i,j+1,k,xsize,ysize)];
				v4 = potential[get_index(i,j-1,k,xsize,ysize)];
				v5 = potential[get_index(i,j,k+1,xsize,ysize)];
				v6 = potential[get_index(i,j,k-1,xsize,ysize)];
				
				currentIndex = get_index(i,j,k,xsize,ysize);
				//printf("%f\n/n",(v1 + v2 + v3 + v4 +v5 +v6 - deltaSquared * source[currentIndex]));
				
				result[currentIndex] = (v1 + v2 + v3 + v4 +v5 +v6 - deltaSquared * source[currentIndex]) / 6;
						
			}	
		}
	}
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
void poisson_dirichlet(double * __restrict__ source,
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
	
	
	
	int startIndex; //Start index for each core...
	int endIndex; // End index for each core...
	int i,n;
	
	double* result = (double *)calloc(xsize * ysize * zsize, sizeof(*potential));


	for (n = 0; n < numiters; n++) 
	{
		std::vector<std::thread> threads;

		// Generates threads equal to the number of cores...
		for (i = 0; i < numcores; ++i)
	    {
			startIndex = (i * xsize / numcores) + 1;
			endIndex = (i + 1) * xsize / numcores + 1;
			threads.push_back(std::thread(poisson_thread_function,startIndex,endIndex,source, potential, xsize, ysize, zsize, delta, result));
		}
		
		
		// Join threads...
		for (auto& thread : threads)
	    {
				thread.join();
		}
		
		// point to new array
		double * temporary_array = result;
		result = potential;
		potential = temporary_array;
		
	}
	
	//free(result);


	
}


