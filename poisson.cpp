#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
#include <thread>
#include <vector>
#include <time.h>
#include <pthread.h>
#include <iostream>

/*
 * poisson.cpp - contains the inital naive_poisson_dirichelt
 * aswell as an optimized function poisson_dirichelt, that utilizes
 * threading, and indexing to maximize caches. 
 * 
 * Authors: Harry Dobbs & Joshua Hulbert
 * 
 * Date: 10th October 2020
 * 
 */

// PLEASE NOTE THE PROVIDED IMPLEMENTATION HAS BEEN ADJUSTED TO USE CALLOC
// THE OUTPUTS OF THE NAIVE AND NEW IMPLEMENTATION HAVE BEEN CROSSED CHECKED
// TO ENSURE EXACT OUTPUTS.


void naive_poisson_dirichlet(double * __restrict__ source,
							double * __restrict__ potential,
							double Vbound,
							unsigned int xsize, unsigned int ysize, unsigned int zsize, double delta,
							unsigned int numiters, unsigned int numcores)
{
	// source[i, j, k] is accessed with source[((k * ysize) + j) * xsize + i]
    // potential[i, j, k] is accessed with potential[((k * ysize) + j) * xsize + i]    
    size_t size = (size_t)ysize * zsize * xsize * sizeof(double);
	double *input = (double *)calloc(size, 1);
	if (!input) {
		fprintf(stderr, "malloc failure\n");
		return;
	}
	for (unsigned int iter = 0; iter < numiters; iter++) {
		for (unsigned int x = 0; x < xsize; x++) {
			for (unsigned int z = 0; z < zsize; z++) { 
				for (unsigned int y = 0; y < ysize; y++) {

					double res = 0;

					if (x < xsize - 1)
						res += input[((z * ysize) + y) * xsize + (x + 1)];
					else
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
		memcpy(input, potential, size);
	}
	free(input);
}


void poisson_thread_function(int start_index, int increment,double * __restrict__ source, double * __restrict__ potential, double Vbound, unsigned int xsize, unsigned int ysize, unsigned int zsize, double delta, double * input)
{
	int x, y, z;
	double res;
	
	for (z = start_index; z <  zsize; z+=increment) {
	  for (y = 0; y < ysize; y++) {
		for (x = 0; x < xsize; x++) {

				res = 0;					
				if (x < xsize - 1){
					res += input[((z * ysize) + y) * xsize + (x + 1)];
					}
				else
					res += Vbound;
					
				if (x > 0){
					res += input[((z * ysize) + y) * xsize + (x - 1)];
					}

				else
					res += Vbound;

				if (y < ysize - 1){
					res += input[((z * ysize) + (y + 1)) * xsize + x];
					}

				else
					res += Vbound;

					
				if (y > 0){
					res += input[((z * ysize) + (y - 1)) * xsize + x];
					}

				else
					res += Vbound;

				if (z < zsize - 1){
					res += input[(((z + 1) * ysize) + y) * xsize + x];
					}
				else
					res += Vbound;
				
				if (z > 0){
					res += input[(((z - 1) * ysize) + y) * xsize + x];					
					}
				else
					res += Vbound;

				res -= delta * delta * source[((z * ysize) + y) * xsize + x];
				res /= 6;
				potential[((z * ysize) + y) * xsize + x] = res;
				
			   }
		   }
	   }
}




/// Solve Poisson's equationinput for a rectangular box with Dirichlet
/// boundary conditions on each face.
/// \param source is a pointer to a flattened 3-D array for the source function
/// \param potential is a pointer to a flattened 3-D array for the calculated potential
/// \param Vbound is the potentipotential = input;al on the boundary
/// \param xsize is the number of elements in the x-direction
/// \param ysize is the number of elements in the y-direction
/// \param zsize is the number of elements in the z-direction
/// \param delta is the voxeinputl spacing in all directions
/// \param numiters is the number of iterations to perform
/// \param numcores is the number of CPU cores to use.  If 0, an optimal number is chosen
void poisson_dirichlet(double * __restrict__ source,double * __restrict__ potential, double Vbound, unsigned int xsize, unsigned int ysize, unsigned int zsize, double delta, unsigned int numiters, unsigned int numcores)
{
	
	size_t size = (size_t)ysize * zsize * xsize * sizeof(double);
	double *input = (double *)calloc(xsize * ysize * zsize, sizeof(*potential));
	double *temp = (double *)calloc(xsize * ysize * zsize, sizeof(*potential));
	if (!input) {
		fprintf(stderr, "malloc failure\n");
		return;
	}
	
	int start_index; //Start index for each core...
	int increment; // End index for each core...
	int i,n;
	
	if(numcores == 0){
		numcores = 1;
	}
	
	std::vector<std::thread> threads;

	
	for (n = 0; n < numiters; n++) 
	{
		
		for(i = 0; i < numcores; i++)
		{	
			start_index = i;
			increment =   numcores ;
			threads.push_back(std::thread(poisson_thread_function,start_index,increment,source, potential, Vbound, xsize, ysize, zsize, delta, input));
		}
		
		for (std::thread & th : threads)
		{
			if (th.joinable())
			{
				th.join();
			} 
		}
	
		
		temp = input;
		input = potential;
		potential = temp;	
	}

	if(numiters % 2 == 0)
	{
		memcpy(potential,input , size);
		free(input);
	}
	else
	{
		free(input);
	}
	
}
	




	


