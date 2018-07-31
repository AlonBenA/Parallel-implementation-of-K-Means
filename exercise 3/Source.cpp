#include<stdio.h>
#include <stdlib.h>
#include<mpi.h>
#include <math.h>
#include<omp.h>

#include "Header.h"

#define FILEWAY "C:/Downloads/Aloninput10.txt" //The file name HAVE to be absolute path + filename and .txt
#define OUTPUTFILEWAY "C:/Downloads/output.txt" //The file name HAVE to be absolute path + filename and .txt
#define MASTER 0
#define ZERO_TAG 0
#define DIDNT_FOUND 0
#define FOUND 1
#define TIME_OUT 2

#define NUMBER_OF_SLAVES 2
#define NUMBER_OF_THREADS 3



struct cluster {
	int id;
	double x;
	double y;
	double diameter;
	int numberOfPoints;
}typedef cluster;

struct Data {
	// number of points
	int NumberOfPoints;
	//number of clusters to find
	int NumberOfClusters;
	// the maximum number of iterations for K - MEAN algorithm.
	int LIMIT;
	//quality measure to stop
	double QM;
	//defines the end of time interval[0, T]
	double T;
	//defines moments t = n*dT, n = { 0, 1, 2, … , T / dT } for which calculate the clusters and the quality
	double dt;

}typedef data;


void initializationOFArray(int array[], int size)
{
	int i;
	for (i = 0; i < size; i++)
	{
		array[i] = 0;
	}
}

float distanceBetweenPointAndCluster(point p1, cluster p2)
{
	//Calculates the distance between a point and cluster
	double  gdistance = ((p2.x - p1.x)*(p2.x - p1.x)) + ((p2.y - p1.y)*(p2.y - p1.y));
	return (float) gdistance;
}

double distanceBetween2Point(point p1, point p2)
{
	//Calculates the distance between 2 points
	double gdistance = ((p2.x - p1.x)*(p2.x - p1.x)) + ((p2.y - p1.y)*(p2.y - p1.y));
	if (gdistance<0)
		return sqrt(-gdistance);
	return sqrt(gdistance);
}

double distanceBetween2Cluster(cluster p1, cluster p2)
{
	//Calculates the distance between 2 cluster
	double gdistance = ((p2.x - p1.x)*(p2.x - p1.x)) + ((p2.y - p1.y)*(p2.y - p1.y));
	if (gdistance<0)
		return sqrt(-gdistance);
	return sqrt(gdistance);
}

void printArrayOfPoints(point *array, int arraySize)
{
	//Print all the informtion of the points
	int i;
	for (i = 0; i < arraySize; i++)
	{
		printf("\npoint number %d Coordinates (%lf,%lf) and vector [%lf,%lf] \n and the closest cluster = %d  \n", i, array[i].x, array[i].y, array[i].Vx, array[i].Vy, array[i].cluster);
	}
}

void printArrayOfClusters(cluster *array, int arraySize)
{
	//Print all the informtion of the Clusters
	int i;
	for (i = 0; i < arraySize; i++)
	{
		printf("cluster number %d, Coordinates(%.5lf, %.5lf), the number of points is %d\n and diameter is %.2lf\n", array[i].id, array[i].x, array[i].y, array[i].numberOfPoints, array[i].diameter);
	}
}

MPI_Datatype createPoint()
{
	//Create the point type for the MPI to recognize 
	MPI_Datatype PointMPIType;
	point p;
	MPI_Datatype type[5] = { MPI_DOUBLE, MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_INT };
	int blocklen[5] = { 1, 1,1,1,1 };
	MPI_Aint disp[5];


	// Create MPI user data type for Point
	disp[0] = (char *)&p.x - (char *)&p;
	disp[1] = (char *)&p.y - (char *)&p;
	disp[2] = (char *)&p.Vx - (char *)&p;
	disp[3] = (char *)&p.Vy - (char *)&p;
	disp[4] = (char *)&p.cluster - (char *)&p;

	MPI_Type_create_struct(5, blocklen, disp, type, &PointMPIType);
	MPI_Type_commit(&PointMPIType);

	return PointMPIType;
}

MPI_Datatype createCluster()
{
	//Create the Cluster type for the MPI to recognize 
	MPI_Datatype clusterMPIType;
	cluster c;
	MPI_Datatype type[5] = { MPI_INT, MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_INT };
	int blocklen[5] = { 1,1,1,1,1 };
	MPI_Aint disp[5];


	// Create MPI user data type for Point
	disp[0] = (char *)&c.id - (char *)&c;
	disp[1] = (char *)&c.x - (char *)&c;
	disp[2] = (char *)&c.y - (char *)&c;
	disp[3] = (char *)&c.diameter - (char *)&c;
	disp[4] = (char *)&c.numberOfPoints - (char *)&c;

	MPI_Type_create_struct(5, blocklen, disp, type, &clusterMPIType);
	MPI_Type_commit(&clusterMPIType);

	return clusterMPIType;
}


MPI_Datatype createData()
{
	//Create the Data type for the MPI to recognize 
	MPI_Datatype DataMPIType;
	data c;
	MPI_Datatype type[6] = { MPI_INT, MPI_INT,MPI_INT,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE };
	int blocklen[6] = { 1,1,1,1,1,1 };
	MPI_Aint disp[6];


	// Create MPI user data type for Data
	disp[0] = (char *)&c.NumberOfPoints - (char *)&c;
	disp[1] = (char *)&c.NumberOfClusters - (char *)&c;
	disp[2] = (char *)&c.LIMIT - (char *)&c;
	disp[3] = (char *)&c.QM - (char *)&c;
	disp[4] = (char *)&c.T - (char *)&c;
	disp[5] = (char *)&c.dt - (char *)&c;

	MPI_Type_create_struct(6, blocklen, disp, type, &DataMPIType);
	MPI_Type_commit(&DataMPIType);

	return DataMPIType;
}

void read_file(point **pointArray, cluster **clusterArray, cluster **tempClusterArray, data *TheData)
{
	FILE *stream;
	char filename[] = FILEWAY;
	int numberOfPoints, k, numberOfCluster, limit;
	double dt, qm, T;
	double x, y, Vx, Vy;


	//Open and get the file handle

	errno_t err = fopen_s(&stream, filename, "r+");
	if (err)
	{
		printf_s("The file was not opened\n");
	}
	else
	{

		// Read data back from file:  
		fscanf_s(stream, "%d %d %lf %lf %d %lf", &numberOfPoints, &numberOfCluster, &T, &dt, &limit, &qm);

		(*TheData).NumberOfPoints = numberOfPoints;
		(*TheData).NumberOfClusters = numberOfCluster;
		(*TheData).T = T;
		(*TheData).dt = dt;
		(*TheData).LIMIT = limit;
		(*TheData).QM = qm;

		//create the arrays for the amount of points and cluster needed
		*pointArray = (point *)malloc(numberOfPoints * sizeof(point));
		*clusterArray = (cluster *)malloc(numberOfCluster * sizeof(cluster));
		*tempClusterArray = (cluster *)malloc(numberOfCluster * sizeof(cluster));

		//read all the points in the file
		for (k = 0; k < numberOfPoints; k++)
		{
			fscanf_s(stream, "%lf%lf%lf%lf", &x, &y, &Vx, &Vy);
			(*pointArray)[k].x = x;
			(*pointArray)[k].y = y;
			(*pointArray)[k].Vx = Vx;
			(*pointArray)[k].Vy = Vy;
			(*pointArray)[k].cluster = -1;

		}

		fclose(stream);
	}
}


void Write_file(cluster *array, int arraySize, double finalTime, double finalQuality)
{
	int i;
	printf("First occurrence at t = %lf with q = %lf", finalTime, finalQuality);
	printf("\n\nend cluster : \n\n");
	printArrayOfClusters(array, arraySize);

	FILE *stream;
	char filename[] = OUTPUTFILEWAY;
	//Open and get the file handle

	errno_t err = fopen_s(&stream, filename, "w+");
	if (err)
	{
		printf_s("The file was not opened\n");
	}
	else
	{
		fprintf(stream, "First occurrence at t = %lf with q = %lf \n", finalTime, finalQuality);
		fprintf(stream, "Centers of the clusters:\n", finalTime, finalQuality);

		//Write all the clusters to the file.
		for (i = 0; i < arraySize; i++)
		{
			fprintf(stream, " %lf %lf \n", array[i].x, array[i].y);
		}

	}

	fclose(stream);

}

void CopyClusterArray(cluster *clusterArray, cluster *tempClusterArray, int NumberOfClusters)
{
	//Copy from one cluster array to the other.
	//if the number of points is zero, the cluster x and y needed to save, so we just update the rest.
	int i;
	for (i = 0; i < NumberOfClusters; i++)
	{
		if (tempClusterArray[i].numberOfPoints != 0)
		{
			clusterArray[i].id = tempClusterArray[i].id;
			clusterArray[i].x = tempClusterArray[i].x;
			clusterArray[i].y = tempClusterArray[i].y;
			clusterArray[i].numberOfPoints = tempClusterArray[i].numberOfPoints;
			clusterArray[i].diameter = tempClusterArray[i].diameter;
		}
		else
		{
			clusterArray[i].id = tempClusterArray[i].id;
			clusterArray[i].numberOfPoints = tempClusterArray[i].numberOfPoints;
			clusterArray[i].diameter = tempClusterArray[i].diameter;
		}
	}
}

void SendToSlave(point *pointArray, cluster *clusterArray, data theData, int NumberOfPoints, int NumberOfClusters, MPI_Datatype PointMPIType, MPI_Datatype clusterMPIType, MPI_Datatype DataMPIType, int NumberOfSlaves)
{
	//Send the points array, cluster Array and all the information needed for the clusters to work separately 
	int i;
	for (i = 1; i <= NumberOfSlaves; i++)
	{
		MPI_Send(&theData, 1, DataMPIType, i, ZERO_TAG, MPI_COMM_WORLD);
		MPI_Send(pointArray, NumberOfPoints, PointMPIType, i, ZERO_TAG, MPI_COMM_WORLD);
		MPI_Send(clusterArray, NumberOfClusters, clusterMPIType, i, ZERO_TAG, MPI_COMM_WORLD);
	}

}

void EndSlaves(int NumberOfSlaves, int SlaveStatus[3], int NumberOfClusters, MPI_Datatype clusterMPIType)
{
	//The End the work of all the threads
	MPI_Status status;
	int EndFlag = FOUND;
	int i;
	int temp;
	int numberOfActiveSlave = NumberOfSlaves;
	double finalData[2] = { 0,0 };
	cluster *tempClusterArray = NULL;
	//Checking how many thread are still active

	for (i = 0; i < NumberOfSlaves; i++)
	{
		if (SlaveStatus[i + 1] != DIDNT_FOUND)
		{
			numberOfActiveSlave -= 1;
		}
	}
	//Gets the message from all the threads that are still active and if they found clusters that are less than the quality
	//the master also gets a message with those clustera and data of who found it and when.
	//but not saveing it, the thread that find the qulity will stop working 
	//and if the thread didn't found qulity the master send him that he found one so he can't stop. 
	for (i = 0; i < numberOfActiveSlave; i++)
	{
		MPI_Recv(&temp, 1, MPI_INT, MPI_ANY_SOURCE, ZERO_TAG, MPI_COMM_WORLD, &status);
		if (temp == FOUND)
		{
			tempClusterArray = (cluster *)malloc(NumberOfClusters * sizeof(cluster));
			MPI_Recv(finalData, 2, MPI_DOUBLE, MPI_ANY_SOURCE, ZERO_TAG, MPI_COMM_WORLD, &status);
			MPI_Recv(tempClusterArray, NumberOfClusters, clusterMPIType, MPI_ANY_SOURCE, ZERO_TAG, MPI_COMM_WORLD, &status);
			SlaveStatus[status.MPI_SOURCE] = FOUND;
		}
		else
		{
			MPI_Send(&EndFlag, 1, MPI_INT, status.MPI_SOURCE, ZERO_TAG, MPI_COMM_WORLD);
		}
	}

}

int CheckWithTheSlaves(int NumberOfSlaves, int SlaveStatus[3], cluster *clusterArray, int NumberOfClusters, MPI_Datatype clusterMPIType, double finalData[2])
{
	MPI_Status status;
	int EndFlag = DIDNT_FOUND;
	int i;
	int temp;
	int numberOfActiveSlave = NumberOfSlaves;
	cluster *tempClusterArray = NULL;

	//Checking how many thread are still active
	for (i = 0; i < NumberOfSlaves; i++)
	{
		if (SlaveStatus[i + 1] != DIDNT_FOUND)
		{
			numberOfActiveSlave -= 1;
		}
	}
	
	//Gets the message from all the threads that are still active and if they are found clusters that are less than the quality
	//the master save the first one and send to the rest that didn't reported the status 
	for (i = 0; i < numberOfActiveSlave; i++)
	{
		MPI_Recv(&temp, 1, MPI_INT, MPI_ANY_SOURCE, ZERO_TAG, MPI_COMM_WORLD, &status);
		if (temp == FOUND)
		{
			tempClusterArray = (cluster *)malloc(NumberOfClusters * sizeof(cluster));
			MPI_Recv(finalData, 2, MPI_DOUBLE, MPI_ANY_SOURCE, ZERO_TAG, MPI_COMM_WORLD, &status);
			MPI_Recv(tempClusterArray, NumberOfClusters, clusterMPIType, MPI_ANY_SOURCE, ZERO_TAG, MPI_COMM_WORLD, &status);
			SlaveStatus[status.MPI_SOURCE] = FOUND;
			if (EndFlag != FOUND)
			{
				EndFlag = FOUND;
				CopyClusterArray(clusterArray, tempClusterArray, NumberOfClusters);
			}
			else
			{

			}
		}
		else
		{
			MPI_Send(&EndFlag, 1, MPI_INT, status.MPI_SOURCE, ZERO_TAG, MPI_COMM_WORLD);
		}
	}
	return EndFlag;
}

void UpdateAllPointsWithCuda(point *pointArray, int NumberOfPoints, double t)
{
	//Update all the points to the right time with cuda
	int workForThread = 100;
	int numberOfThreads = (NumberOfPoints / workForThread);
	int numberOfThreadsForBlock;
	int numberOfBlocks;

	if (numberOfThreads > 1000)
	{
		numberOfBlocks = numberOfThreads / 1000;
		numberOfThreadsForBlock = 1000;
		if (numberOfThreads % 1000 != 0)
		{
			numberOfBlocks++;
		}
	}
	else
	{
		//if there is less then 100 points there is only 1 thread in cuda.
		if (numberOfThreads != 0)
		{
			numberOfThreadsForBlock = numberOfThreads;
		}
		else
		{
			numberOfThreads = 1;
			numberOfThreadsForBlock = 1;
		}

		numberOfBlocks = 1;
	}


	cudaError_t cudaStatus = CudaMain(pointArray,NumberOfPoints,t, numberOfBlocks,numberOfThreadsForBlock,workForThread);

	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "addWithCuda failed!");
	}
	
}

void SelectFirstClusters(point *pointArray, cluster *clusterArray, int NumberOfPoints, int NumberOfClusters)
{
	int i;
	for (i = 0; i< NumberOfClusters; i++)
	{
		clusterArray[i].id = i;
		clusterArray[i].x = pointArray[i].x;
		clusterArray[i].y = pointArray[i].y;
		clusterArray[i].numberOfPoints = 0;
		clusterArray[i].diameter = 0;
	}
}

int ClosetClustersForEachPoint(point *pointArray, cluster *clusterArray, int NumberOfPoints, int NumberOfClusters)
{
	//for each point Calculates the nearest cluster, if it is different from what it used to be
	//then breakFlag is Changing and another round is needed. 
	int i, j, tempId;
	float min, temp;
	int breakFlag = 0;


	for (i = 0; i < NumberOfPoints; i++)
	{
		min = distanceBetweenPointAndCluster(pointArray[i], clusterArray[0]);
		tempId = 0;

		for (j = 1; j < NumberOfClusters; j++)
		{
			temp = distanceBetweenPointAndCluster(pointArray[i], clusterArray[j]);

			if (temp < min)
			{
				tempId = j;
				min = temp;
			}
		}

		if (pointArray[i].cluster != tempId)
		{
			breakFlag = 1;
		}
		pointArray[i].cluster = tempId;
	}

	return breakFlag;
}

void createNewClustersArray(point *pointArray, cluster *NewClusterArray, int NumberOfPoints, int NumberOfClusters)
{
	int i, j, clusterId;
	double newX, newY, NumberOfPointInCluster;

	for (i = 0; i < NumberOfClusters; i++)
	{
		NewClusterArray[i].numberOfPoints = 0;
		NewClusterArray[i].x = 0;
		NewClusterArray[i].y = 0;
		NewClusterArray[i].id = i;
		NewClusterArray[i].diameter = 0;
	}
	//omp
#pragma omp parallel for shared(pointArray,NumberOfPoints) private(clusterId,i)
	for (i = 0; i < NumberOfPoints; i++)
	{
		clusterId = pointArray[i].cluster;
		NewClusterArray[clusterId].x += pointArray[i].x;
		NewClusterArray[clusterId].y += pointArray[i].y;
		NewClusterArray[clusterId].numberOfPoints++;
	}

	for (i = 0; i < NumberOfClusters; i++)
	{
		if (NewClusterArray[i].numberOfPoints != 0)
		{
			NewClusterArray[i].x = NewClusterArray[i].x / NewClusterArray[i].numberOfPoints;
			NewClusterArray[i].y = NewClusterArray[i].y / NewClusterArray[i].numberOfPoints;
		}
	}
}

void diameterForClusters(point *pointArray, cluster *ClusterArray, int NumberOfPoints, int NumberOfClusters)
{
	int i, j, k;
	double diameter;
	//omp
#pragma omp parallel for shared(pointArray,ClusterArray,NumberOfPoints,NumberOfClusters) private(j,k,diameter,i)
	for (i = 0; i < NumberOfClusters; i++)
	{
		ClusterArray[i].diameter = 0;
		diameter = 0;
		for (j = 0; j < NumberOfPoints - 1; j++)
		{
			if (pointArray[j].cluster == i)
			{
				for (k = j + 1; k < NumberOfPoints; k++)
				{
					if (pointArray[k].cluster == i)
					{
						diameter = distanceBetween2Point(pointArray[j], pointArray[k]);
						if (ClusterArray[i].diameter < diameter)
							ClusterArray[i].diameter = diameter;
					}
				}
			}

		}
	}


}

double qualityCalculate(cluster *ClusterArray, int NumberOfClusters)
{
	double quality = 0;
	double distance;
	int i, j;
	int div = (NumberOfClusters*(NumberOfClusters - 1));
	for (i = 0; i < NumberOfClusters; i++)
	{
		for (j = 0; j < NumberOfClusters; j++)
		{
			if (i != j)
			{
				distance = distanceBetween2Cluster(ClusterArray[i], ClusterArray[j]);
				quality += ClusterArray[i].diameter / distance;
			}
		}

	}

	quality = quality / div;
	return quality;
}

double Kmean(point *pointArray, cluster *clusterArray, cluster *tempClusterArray, int NumberOfPoints, int NumberOfClusters, int limit)
{
	int i = 0;
	int breakFlag;
	double quality;

	//Select K CLUSTERS(the first k points)
	SelectFirstClusters(pointArray, clusterArray, NumberOfPoints, NumberOfClusters);


	breakFlag = 1;
	i = 0;

	while (breakFlag && i < limit)
	{
		if (i == 0)
		{
			//for each point find the closest CLUSTER
			ClosetClustersForEachPoint(pointArray, clusterArray, NumberOfPoints, NumberOfClusters);
		}
		else
		{
			breakFlag = ClosetClustersForEachPoint(pointArray, tempClusterArray, NumberOfPoints, NumberOfClusters);
		}


		createNewClustersArray(pointArray, tempClusterArray, NumberOfPoints, NumberOfClusters);

		CopyClusterArray(clusterArray, tempClusterArray, NumberOfClusters);
		i++;
	}


	diameterForClusters(pointArray, clusterArray, NumberOfPoints, NumberOfClusters);
	quality = qualityCalculate(clusterArray, NumberOfClusters);

	return quality;
}



void MasterWork(MPI_Datatype PointMPIType, MPI_Datatype clusterMPIType, MPI_Datatype DataMPIType, int myid)
{
	double t = 0;
	double finalTime, finalQuality;
	double finalData[2] = { 0,0 };
	int SlaveStatus[NUMBER_OF_THREADS];
	int Endflag = 0;
	double quality;
	int NumberOfPoints;
	int NumberOfClusters;
	point *pointArray = NULL;
	cluster *clusterArray = NULL;
	cluster *tempClusterArray = NULL;
	data theData;
	MPI_Status status;


	initializationOFArray(SlaveStatus, NUMBER_OF_THREADS);
	read_file(&pointArray, &clusterArray, &tempClusterArray, &theData);

	NumberOfPoints = theData.NumberOfPoints;
	NumberOfClusters = theData.NumberOfClusters;


	SendToSlave(pointArray, clusterArray, theData, NumberOfPoints, NumberOfClusters, PointMPIType, clusterMPIType, DataMPIType, NUMBER_OF_SLAVES);

	for (t = myid*theData.dt; t < theData.T; t += NUMBER_OF_THREADS*theData.dt)
	{

		if (t == 0)
		{

		}
		else
		{

			UpdateAllPointsWithCuda(pointArray, NumberOfPoints, NUMBER_OF_THREADS*theData.dt);

		}

	

		//kmean
		quality = Kmean(pointArray, clusterArray, tempClusterArray, NumberOfPoints, NumberOfClusters, theData.LIMIT);



		if (quality < theData.QM)
		{

			finalTime = t;
			finalQuality = quality;

			EndSlaves(NUMBER_OF_SLAVES, SlaveStatus, NumberOfClusters, clusterMPIType);
			break;
		}
		else
		{
			Endflag = CheckWithTheSlaves(NUMBER_OF_SLAVES, SlaveStatus, clusterArray, NumberOfClusters, clusterMPIType, finalData);
		}

		if (Endflag == FOUND)
		{

			finalTime = finalData[0];
			finalQuality = finalData[1];
			break;
		}

	}

	if (t >= theData.T)
	{
		finalTime = -myid;
		finalQuality = -myid;
		EndSlaves(NUMBER_OF_SLAVES, SlaveStatus, NumberOfClusters, clusterMPIType);
	}

	if (finalQuality != 0)
	{
		Write_file(clusterArray, NumberOfClusters, finalTime, finalQuality);
	}
	else
	{
		Write_file(clusterArray, NumberOfClusters, t, quality);
		printf("\n\n never occurrence \n\n");
	}
}

void SlaveWork(MPI_Datatype PointMPIType, MPI_Datatype clusterMPIType, MPI_Datatype DataMPIType, int myid)
{
	int EndFlag = DIDNT_FOUND;
	double t = 0;
	double quality;
	double finalData[2] = { 0,0 };
	int NumberOfPoints, NumberOfClusters;
	MPI_Status status;
	point *pointArray = NULL;
	cluster *clusterArray = NULL;
	cluster *tempClusterArray = NULL;
	data theData;

	//Getting all the information from the master
	MPI_Recv(&theData, 1, DataMPIType, MASTER, ZERO_TAG, MPI_COMM_WORLD, &status);
	NumberOfPoints = theData.NumberOfPoints;
	NumberOfClusters = theData.NumberOfClusters;

	pointArray = (point *)malloc(NumberOfPoints * sizeof(point));
	clusterArray = (cluster *)malloc(NumberOfClusters * sizeof(cluster));
	tempClusterArray = (cluster *)malloc(NumberOfClusters * sizeof(cluster));

	MPI_Recv(pointArray, NumberOfPoints, PointMPIType, MASTER, ZERO_TAG, MPI_COMM_WORLD, &status);
	MPI_Recv(clusterArray, NumberOfClusters, clusterMPIType, MASTER, ZERO_TAG, MPI_COMM_WORLD, &status);

	for (t = myid*theData.dt; t < theData.T; t += NUMBER_OF_THREADS*theData.dt)
	{
		if (t == myid*theData.dt)
		{
			
			UpdateAllPointsWithCuda(pointArray, NumberOfPoints, t);
		
		}
		else
		{
			
			UpdateAllPointsWithCuda(pointArray, NumberOfPoints, NUMBER_OF_THREADS*theData.dt);
			
		}

		//kmean
		quality = Kmean(pointArray, clusterArray, tempClusterArray, NumberOfPoints, NumberOfClusters, theData.LIMIT);


		if (quality < theData.QM)
		{

			EndFlag = FOUND;
			MPI_Send(&EndFlag, 1, MPI_INT, MASTER, ZERO_TAG, MPI_COMM_WORLD);

			finalData[0] = t;
			finalData[1] = quality;
			MPI_Send(finalData, 2, MPI_DOUBLE, MASTER, ZERO_TAG, MPI_COMM_WORLD);

			MPI_Send(clusterArray, NumberOfClusters, clusterMPIType, MASTER, ZERO_TAG, MPI_COMM_WORLD);
			break;
		}
		else
		{
			EndFlag = DIDNT_FOUND;
			MPI_Send(&EndFlag, 1, MPI_INT, MASTER, ZERO_TAG, MPI_COMM_WORLD);
			MPI_Recv(&EndFlag, 1, MPI_INT, MASTER, ZERO_TAG, MPI_COMM_WORLD, &status);
		}

		if (EndFlag == FOUND)
		{
			break;
		}

	}

	printf("my id = %d t = %lf \n", myid, t);

	if (t >= theData.T)
	{
		EndFlag = TIME_OUT;
		MPI_Send(&EndFlag, 1, MPI_INT, MASTER, ZERO_TAG, MPI_COMM_WORLD);
	}

}

int main(int argc, char *argv[])
{
	int  namelen, numprocs, myid;
	char processor_name[MPI_MAX_PROCESSOR_NAME];
	MPI_Datatype PointMPIType;
	MPI_Datatype clusterMPIType;
	MPI_Datatype DataMPIType;
	MPI_Init(&argc, &argv);
	MPI_Status status;
	double t1, t2;

	MPI_Comm_rank(MPI_COMM_WORLD, &myid);

	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

	MPI_Get_processor_name(processor_name, &namelen);

	if (myid == 0)
		t1 = MPI_Wtime();

	PointMPIType = createPoint();
	clusterMPIType = createCluster();
	DataMPIType = createData();


	if (myid == 0) {
		MasterWork(PointMPIType, clusterMPIType, DataMPIType, myid);
	}
	else {
		SlaveWork(PointMPIType, clusterMPIType, DataMPIType, myid);

	}

	if (myid == 0)
	{
		t2 = MPI_Wtime();
		printf("MPI_Wtime measured to be: %3.3f\n", t2 - t1);
		fflush(stdout);
	}

	MPI_Finalize();
	return 0;
}
