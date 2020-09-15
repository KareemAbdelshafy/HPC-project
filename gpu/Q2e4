#include <cuda.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define BLOCK_SIZE 1024 // You can change this
#define NUM_OF_ELEMS 4096 // You can change this


__global__  void total(float * input, float * output, int len) 
{
    int tid = blockIdx.x*blockDim.x + threadIdx.x;

    // Traverse reduction tree
   for (unsigned int stride = NUM_OF_ELEMS/2; stride > 0; stride /= 2)
    {
        if ((tid < stride) && (input[tid] > input[tid + stride]) ){
           input[tid] = input[tid + stride];}

    __syncthreads();
    }   

    // Write the computed sum of the block to the output vector at correct index
    if (tid == 0)
    {
        output[0] = input[0];
    }
}

int main(int argc, char ** argv) 
{
    float * hostInput; // The input 1D vector
    float * hostOutput; // The output vector
    float * deviceInput;
    float * deviceOutput;

    int numInputElements = NUM_OF_ELEMS; // number of elements in the input list
    //int numOutputElements; // number of elements in the output list
    hostInput = (float *) malloc(sizeof(float) * numInputElements);

    srand(time(NULL));
    for (int i=0; i < NUM_OF_ELEMS; i++)
    {
        hostInput[i] = rand();  
    }

    printf("host %f and %f and %f \n", hostInput[10] , hostInput[20] , hostInput[30]);


    hostOutput = (float*) malloc(sizeof(float));

    //@@ Allocate GPU memory here
    cudaMalloc((void **)&deviceInput, numInputElements * sizeof(float));
    cudaMalloc((void **)&deviceOutput, sizeof(float));

    // Copy memory to the GPU here
    cudaMemcpy(deviceInput, hostInput, numInputElements * sizeof(float), cudaMemcpyHostToDevice);

    // Initialize the grid and block dimensions here
    dim3 DimGrid( 1, 1, 1);
    dim3 DimBlock(BLOCK_SIZE, 1, 1);

    // Launch the GPU Kernel here
    total<<<DimGrid, DimBlock>>>(deviceInput, deviceOutput, numInputElements);

    // Copy the GPU memory back to the CPU here
    cudaMemcpy(hostOutput, deviceOutput, sizeof(float), cudaMemcpyDeviceToHost);

     printf("Reduced Sum from GPU = %f \n", hostOutput[0]);  
     //printf("Reduced Sum from GPU = %d last value %f \n", numOutputElements);   

    // Free the GPU memory here
    cudaFree(deviceInput);
    cudaFree(deviceOutput); 
    free(hostInput);
    free(hostOutput);

    return 0;
}

