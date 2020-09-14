#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <cuda.h>
#include <math.h>

__host__ void init_vects(int vect_len,float *h_vect1,float *h_vect2);
__global__ void vec_add(int vect_len, float *d_vect1, float *d_vect2, float *d_sum);

int main(int argc,char **argv)
{
cudaEvent_t start=0;
cudaEvent_t stop=0;
float time=0;
cudaEventCreate(&start);
cudaEventCreate(&stop);

int vect_len=pow(2,18);
float h_vect1[vect_len], h_vect2[vect_len], h_sum[vect_len];
float *d_vect1, *d_vect2, *d_sum;
// initialization
init_vects(vect_len, h_vect1, h_vect2);

// tranfer vectors to global memory
cudaMalloc((void **)&d_vect1 , vect_len*sizeof(float) );
cudaMalloc((void **)&d_vect2 , vect_len*sizeof(float) );
cudaMalloc((void **)&d_sum , vect_len*sizeof(float) );

cudaMemcpy (d_vect1 , h_vect1 , vect_len*sizeof(float) , cudaMemcpyHostToDevice);
cudaMemcpy (d_vect2 , h_vect2 , vect_len*sizeof(float) , cudaMemcpyHostToDevice);

// determine block and grid size.

dim3 DimGrid((vect_len/1024),1 ,1);
dim3 DimBlock(1024,1,1);


cudaEventRecord(start,0);
vec_add<<<DimGrid,DimBlock>>>(vect_len, d_vect1 ,d_vect2 ,d_sum);
cudaEventRecord(stop,0);
cudaEventSynchronize(stop);

cudaMemcpy(h_sum , d_sum , vect_len*sizeof(float) , cudaMemcpyDeviceToHost);
//Free the Device array
cudaFree (d_vect1);
cudaFree (d_vect2);
cudaFree (d_sum);
cudaEventElapsedTime(&time,start,stop);
printf("time of the Kernel %f \n",time );

printf("v1=%f ,, v2 =%f ,, sum=%f \n",h_vect1[0],h_vect2[0],h_sum[0]);

return 0;
}

__global__ void vec_add(int vect_len, float *d_vect1, float *d_vect2, float *d_sum){

int tid = blockIdx.x * blockDim.x + threadIdx.x ;

//if(tid<vect_len)
d_sum[tid]= d_vect1[tid] + d_vect2[tid];

}


__host__ void init_vects(int vect_len,float *h_vect1,float *h_vect2){
srand(time(NULL));
 for (int i=0; i<vect_len; i++){
    h_vect1[i] = rand();
    h_vect2[i] = rand();
  }
}
