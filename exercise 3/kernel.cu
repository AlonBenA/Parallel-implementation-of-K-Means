
#include <stdio.h>
#include "Header.h"


__global__ void updatePointsKernel(point *pointArray,int NumberOfPoints ,const int workForThread, const int numberOfThreadsForBlock,const double t)
{
	int myId = threadIdx.x + blockIdx.x * numberOfThreadsForBlock;
	int start = myId * workForThread;
	int end = start + workForThread;
	int i;
	if(end > NumberOfPoints )
		end = NumberOfPoints ;

	for (i = start; i < end; i++)
	{
		pointArray[i].x += pointArray[i].Vx*t;
		pointArray[i].y += pointArray[i].Vy*t;
	}

}


// Helper function for using CUDA to update points in parallel.
cudaError_t CudaMain(point *pointArray, int NumberOfPoints,double t, int numberOfBlocks,int numberOfThreadsForBlock, int workForThread)
{
point *dev_pointsArray = 0;
    cudaError_t cudaStatus;

    // Choose which GPU to run on, change this on a multi-GPU system.
    cudaStatus = cudaSetDevice(0);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");

    }

    // Allocate GPU buffers for the points (two input, one output)    .
    cudaStatus = cudaMalloc((void**)&dev_pointsArray, NumberOfPoints * sizeof(point));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
    }

    // Copy the points from host memory to GPU buffers.
    cudaStatus = cudaMemcpy(dev_pointsArray, pointArray, NumberOfPoints * sizeof(point), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");

    }


    // Launch a kernel on the GPU 
    updatePointsKernel<<<numberOfBlocks, numberOfThreadsForBlock >>>(dev_pointsArray,NumberOfPoints ,workForThread, numberOfThreadsForBlock, t);

    // Check for any errors launching the kernel
    cudaStatus = cudaGetLastError();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "addKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));

    }
    
    // cudaDeviceSynchronize waits for the kernel to finish, and returns
    // any errors encountered during the launch.
    cudaStatus = cudaDeviceSynchronize();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching addKernel!\n", cudaStatus);

    }

    // Copy the update points from GPU buffer to host memory.
    cudaStatus = cudaMemcpy(pointArray, dev_pointsArray, NumberOfPoints * sizeof(point), cudaMemcpyDeviceToHost);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");

    }


    cudaFree(dev_pointsArray);

    
    return cudaStatus;
}
