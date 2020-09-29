#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
#include <thread>
#include <vector>
#include <time.h>
#include <pthread.h>



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


void poisson_thread_function(int start_index, int end_index,double * __restrict__ source, double * __restrict__ potential, double Vbound, unsigned int xsize, unsigned int ysize, unsigned int zsize, double delta, double * input)
{
	unsigned int x, y, z;
	double res;
	

	   
	  for (z = 0; z < zsize; z++) {
		for (y = 0; y < ysize; y++) {
			for (x = start_index ; x < end_index; x+= xsize-1) {
				
				res = 0;
				
				if (x < xsize - 1)
					res += input[((z * ysize) + y) * xsize + (x + 1)];
				else
					printf("entered last one");
					res += Vbound;
					
				if (x > 0)
					res += input[((z * ysize) + y) * xsize + (x - 1)];
				else
					res += Vbound;

				if (y < ysize - 1)
					res += input[((z * ysize) + (y + 1)) * xsize + x];
				else
					res += Vbound;
					
				if (y > 0)
					res += input[((z * ysize) + (y - 1)) * xsize + x];
				else
					res += Vbound;

				if (z < zsize - 1)
					res += input[(((z + 1) * ysize) + y) * xsize + x];
				else
					res += Vbound;
				
				if (z > 0)
					res += input[(((z - 1) * ysize) + y) * xsize + x];
				else
					res += Vbound;

				res -= delta * delta * source[((z * ysize) + y) * xsize + x];
				res /= 6;
				potential[((z * ysize) + y) * xsize + x] = res;
				
			   }
		   }
	   }
	   
	   
	   
	   for (z = 1; z < zsize-1; z++) {
		for (y = 1; y < ysize-1; y++) {
			for (x = start_index+1; x < end_index-1; x++) {
				//printf("entered");
				res = 0;
				res += input[((z * ysize) + y) * xsize + (x + 1)];
				res += input[((z * ysize) + y) * xsize + (x - 1)];
				res += input[((z * ysize) + (y + 1)) * xsize + x];
				res += input[((z * ysize) + (y - 1)) * xsize + x];
				res += input[(((z + 1) * ysize) + y) * xsize + x];
				res += input[(((z - 1) * ysize) + y) * xsize + x];
				res -= delta * delta * source[((z * ysize) + y) * xsize + x];
				res /= 6;
				potential[((z * ysize) + y) * xsize + x] = res;
				
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
void poisson_dirichlet(double * __restrict__ source,double * __restrict__ potential, double Vbound, unsigned int xsize, unsigned int ysize, unsigned int zsize, double delta, unsigned int numiters, unsigned int numcores)
{
	
	size_t size = (size_t)ysize * zsize * xsize * sizeof(double);
	double *input = (double *)malloc(size);
	if (!input) {
		fprintf(stderr, "malloc failure\n");
		return;
	}
	
	int start_index; //Start index for each core...
	int end_index; // End index for each core...
	int i,n;
	
	memcpy(input, source, size);

	if(numcores == 0){
		numcores = 1;
	}
	
	std::vector<std::thread> threads;

	
	for (n = 0; n < numiters; n++) 
	{
		
		for(i = 0; i < numcores; i++)
		{	

			start_index = (i * xsize / numcores);
			end_index = (i + 1) * xsize / numcores ;
			threads.push_back(std::thread(poisson_thread_function,start_index,end_index,source, potential, 0, xsize, ysize, zsize, delta, input));
		
		}
		
		for (std::thread & th : threads)
		{
			if (th.joinable())
				th.join();
		}
		
		memcpy(input, potential, size);
	}
	
	free(input);
}
	




	


