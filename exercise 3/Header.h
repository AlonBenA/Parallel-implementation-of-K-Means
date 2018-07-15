#pragma once
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

struct Point {
	double x;
	double y;
	double Vx;
	double Vy;
	int cluster;
}typedef point;

cudaError_t CudaMain(point *pointArray, int NumberOfPoints, double t, int numberOfBlocks, int numberOfThreadsForBlock, int workForThread);