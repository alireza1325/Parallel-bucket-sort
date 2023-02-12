// lab3.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include "Utilities.h"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <random>

using namespace std;


template<typename T>
T random(T range_from, T range_to) {
	std::random_device                  rand_dev;
	std::mt19937                        generator(rand_dev());
	std::uniform_int_distribution<T>    distr(range_from, range_to);
	return distr(generator);
}


void parallelRange(int globalstart, int globalstop, int irank, int nproc, int& localstart, int& localstop, int& localcount)
{
	int nrows = globalstop - globalstart + 1;
	int divisor = nrows / nproc;
	int remainder = nrows % nproc;
	int offset;
	if (irank < remainder) offset = irank;
	else offset = remainder;

	localstart = irank * divisor + globalstart + offset;
	localstop = localstart + divisor - 1;
	if (remainder > irank) localstop += 1;
	localcount = localstop - localstart + 1;
}



template<class T>
std::vector<T> bucketSort(T globalMin, T globalMax, int nproc, vector<T>localNums)
{
	// divide the recieved data into nproc buckets 
	int bucketRank;
	vector<vector<T>> sendBuckets = vector<vector<T>>(nproc);
	for (int i = 0; i < localNums.size(); i++)
	{
		bucketRank = floor(nproc * (localNums[i] - globalMin) / (globalMax - globalMin)); // normalize the data into range [0, nproc]
		if (bucketRank == nproc)
		{
			bucketRank = nproc - 1;
		}
		sendBuckets[bucketRank].push_back(localNums[i]);
	}
	//  Distribute buckest to all processors and receive my bucket to sort
	vector<vector<T>> recvBuckets = vector<vector<T>>();
	MPI_Alltoall_vecvecT(sendBuckets, recvBuckets);

	// convert myBucket which is a vector of vector to a long one dimentioanl array (vector)
	vector<int> myBucket = vector<int>();
	for (int bucket = 0; bucket < recvBuckets.size(); bucket++) {
		for (int i = 0; i < recvBuckets[bucket].size(); i++) {
			myBucket.push_back(recvBuckets[bucket][i]);
		}
	}
	// Sort the bucket sequentially
	std::sort(myBucket.begin(), myBucket.end());

	// Return the sorted bucket 
	return myBucket;
}


template<class T>
T findMedian(int m, vector<T>sortedBucket) {
	int rank = 0;
	int nproc = 1;

	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	// gather all bucketsizes
	int myLength = sortedBucket.size();
	std::vector<int> allLengths(nproc);
	MPI_Gather(&myLength, 1, MPI_INT, &allLengths[0], 1, MPI_INT, 0, MPI_COMM_WORLD);

	// find the proceesor that has the desired median (m) and its corresponding index in that processor 
	int buff = 0;
	int winRank = -1;
	int localM;

	if (rank == 0) {
		for (int i = 0; i < nproc; i++) {
			if (m <= buff + allLengths[i]) {
				winRank = i;
				localM = m - buff-1;
				break;
			}
			else buff += allLengths[i];
		}
		assert(winRank >= 0); // if the entered number m is out of range, the program will be aborted
	}

	// broadcast the winner rank to all processor
	MPI_Bcast(&winRank, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
	// broadcast the local M to all processor
	MPI_Bcast(&localM, 1, MPI_INT, 0, MPI_COMM_WORLD);

	T median;
	if (winRank == rank) {
		median = sortedBucket[localM];
		std::cout << "I am rank " << rank << " and I got the desired median which is " << median << std::endl;
	}
	
	// let al processor knows what is median value
	MPI_Bcast(&median, sizeof(T), MPI_BYTE, winRank, MPI_COMM_WORLD);
	return median;
}

template<class T>
void data_to_file(string file_name,std::vector<T> data)
{
	ofstream out_to_file;
	// clear the file content if it is already existed 
	out_to_file.open(file_name, std::ofstream::out | std::ofstream::trunc);
	out_to_file.close();
	// Store the data
	out_to_file.open(file_name, std::ofstream::out | std::ofstream::app);
	for (int i = 0; i < data.size(); i++) out_to_file << data[i] << std::endl;
	out_to_file.close();
}

int main(int argc, char** argv)
{

	int rank, nproc;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (argc < 2)
	{
		std::cout << "You need to define median" << std::endl;
	}

	int m = stoi(argv[1]);

	if (rank == 0)
	{
		
		int N = 1000;  // number of random numbers
		int minVal = 0;  // minimum value of random numbers	
		int maxVal = 99;	// maximum value of random numbers
		std::vector<int> randNums;
		//randNums = random(0, 84);
		srand(time(0));
		// generate N uniformly distributed random numbers in region of [0,maxVal]
		for (int i = 0; i < N; i++) {
			//double randNum =  random(minVal, maxVal);
			double randNum =  (double)rand() / (double)RAND_MAX;
			randNum = (maxVal - minVal) * randNum + minVal;
			randNums.push_back(randNum);
		}
		// display them
		//std::cout << "Generated Random numbers in processor zero: "<<std::endl;
		//for (int i = 0; i < randNums.size(); i++) std::cout << randNums[i] << ",\t";
		std::cout << std::endl;
		int MaxNum = *max_element(randNums.begin(), randNums.end()); // find the global maximum 
		int MinNum = *min_element(randNums.begin(), randNums.end()); // find the global minimum
		std::cout << "Maximum number " << MaxNum <<std::endl;
		std::cout << "Minimum number " << MinNum << std::endl;

		
		// find the amount of random numbers which is going to be sent to each processor
		int datainfo[3] = { MinNum,MaxNum,0 };
		
		std::vector<int> localstarts(nproc), localstops(nproc), localcounts(nproc), recvbuf;
		for (int irank = 0; irank < nproc; irank++)
		{
			parallelRange(0, N-1, irank, nproc, localstarts[irank], localstops[irank], localcounts[irank]);
			int sendcounts = localcounts[irank];
			if (irank != 0)
			{
				datainfo[2] = sendcounts;
				// send data information to other processors
				MPI_Send(&datainfo, 3, MPI_INT, irank, 99, MPI_COMM_WORLD);
			}
		}
		int recvcounts = localcounts[0];
		recvbuf.resize(recvcounts);
		
		// distribute each processor share ( in this section I didn't create buckets. I just send roughly same number of random numbers to each processor
		MPI_Scatterv(&randNums[0], &localcounts[0], &localstarts[0], MPI_INT, &recvbuf[0], recvcounts, MPI_INT, 0, MPI_COMM_WORLD);

		// start timing 
		MPI_Barrier(MPI_COMM_WORLD);
		double start_time = MPI_Wtime();

		// Do the bucket sort 
		vector<int> mySortedBucket = vector<int>();
		mySortedBucket = bucketSort(MinNum, MaxNum, nproc, recvbuf);

		// finish timing 
		MPI_Barrier(MPI_COMM_WORLD);
		double end_time = MPI_Wtime();

		// Find the median value and processor that has it
		int desiredMedian = findMedian(m, mySortedBucket);


		// Store received and sorted numbers into txt files
		stringstream ss;
		ss << rank;
		string fname = "Received numbers in processor " + ss.str() + " .txt";
		data_to_file(fname, recvbuf);
		
		fname = "Sorted numbers in processor " + ss.str() + " .txt";
		data_to_file(fname, mySortedBucket);

		std::cout<< "Parallel time = "<< end_time - start_time << " seconds." << std::endl;
	}
	else
	{
		

		int datainfo[3];
		// receive required information
		MPI_Recv(&datainfo, 3, MPI_INT, 0, 99, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		
		int MinNum = datainfo[0];
		int MaxNum = datainfo[1];
		int recvcounts = datainfo[2];

		std::vector<int> recvbuf(recvcounts);
		MPI_Barrier(MPI_COMM_WORLD);

		MPI_Scatterv(NULL,NULL, NULL,NULL, &recvbuf[0], recvcounts, MPI_INT, 0, MPI_COMM_WORLD);
		//std::cout << "Received numbers in processor " << rank << std::endl;
		
		vector<int> mySortedBucket = vector<int>();
		mySortedBucket = bucketSort(MinNum, MaxNum, nproc, recvbuf);

		MPI_Barrier(MPI_COMM_WORLD);
		int desiredMedian = findMedian(m, mySortedBucket);



		// Store received and sorted numbers into txt files
		stringstream ss;
		ss << rank;
		string fname = "Received numbers in processor " + ss.str() + " .txt";
		data_to_file(fname, recvbuf);

		fname = "Sorted numbers in processor " + ss.str() + " .txt";
		data_to_file(fname, mySortedBucket);
	}

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
}
