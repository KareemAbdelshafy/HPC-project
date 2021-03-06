#include <cuda.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define BLOCK_SIZE 1024 // You can change this
//#define NUM_OF_ELEMS 1e6 // You can change this


__global__  void total(float * input, float * output, int len) 
{
    int tid = blockIdx.x*blockDim.x + threadIdx.x;

    // Traverse reduction tree
   for (unsigned int stride = len/2; stride > 0; stride /= 2)
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
    int NUM_OF_ELEMS=pow(2,18);

cudaEvent_t start=0;
cudaEvent_t stop=0;
float timef=0;
cudaEventCreate(&start);
cudaEventCreate(&stop);

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
  

    dim3 DimGrid(ceil(NUM_OF_ELEMS/1024), 1, 1);
    dim3 DimBlock(BLOCK_SIZE, 1, 1);

    // Launch the GPU Kernel here
    cudaEventRecord(start,0);
    total<<<DimGrid, DimBlock>>>(deviceInput, deviceOutput, numInputElements);
    cudaEventRecord(stop,0);
    cudaEventSynchronize(stop);

    // Copy the GPU memory back to the CPU here
    cudaMemcpy(hostOutput, deviceOutput, sizeof(float), cudaMemcpyDeviceToHost);

     printf("Reduced Sum from GPU = %f \n", hostOutput[0]);  
     cudaEventElapsedTime(&timef,start,stop);
     printf("time of the Kernel %f \n",timef );  

    // Free the GPU memory here
    cudaFree(deviceInput);
    cudaFree(deviceOutput); 
    free(hostInput);
    free(hostOutput);

    return 0;
}

