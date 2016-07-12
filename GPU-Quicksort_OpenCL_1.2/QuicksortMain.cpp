/* ************************************************************************* *\
               INTEL CORPORATION PROPRIETARY INFORMATION
     This software is supplied under the terms of a license agreement or 
     nondisclosure agreement with Intel Corporation and may not be copied 
     or disclosed except in accordance with the terms of that agreement. 
        Copyright (C) 2014 Intel Corporation. All Rights Reserved.
\* ************************************************************************* */

// QuicksortMain.cpp : Defines the entry point for the console application.
//

#include <stdio.h>
#include <windows.h>
#include <tchar.h>
#include <assert.h>
#include <string.h>

#include "OpenCLUtils.h"
#include <math.h>
#include <iostream>
#include <algorithm>
#include <iterator>
#include <vector>
#include <map>

// Types:
typedef unsigned int uint;
typedef	uint QuicksortFlag;

#define READ_ALIGNMENT  4096 // Intel recommended alignment
#define WRITE_ALIGNMENT 4096 // Intel recommended alignment

#define MAX_TEST_KERNELS 100
typedef struct
{
	char				pKernelName[256];
	uint				numBytesPerWorkItemRead;
	cl_kernel			kernelHdl;
	QuicksortFlag			testFlag;
} QuicksortKernel;

typedef struct
{	
	// CL platform handles:
	cl_device_id		deviceID;
	cl_context			contextHdl;
	cl_program			programHdl;
	cl_command_queue	cmdQHdl;

	cl_mem buffIn;
	cl_mem buffOut;
	
	QuicksortKernel	    	pTests[MAX_TEST_KERNELS];
	uint				numTests;
} OCLResources;

// Globals:
cl_int		ciErrNum;

static cl_kernel gqsort_kernel, lqsort_kernel;

void Cleanup(OCLResources* pOCL, int iExitCode, bool bExit, char* optionalErrorMessage)
{
	if (optionalErrorMessage) printf ("%s\n", optionalErrorMessage);

	memset(pOCL, 0, sizeof (OCLResources));

	if (pOCL->programHdl)		{ clReleaseProgram(pOCL->programHdl);		pOCL->programHdl=NULL;	}
	if (pOCL->cmdQHdl)			{ clReleaseCommandQueue(pOCL->cmdQHdl);		pOCL->cmdQHdl=NULL;		}
	if (pOCL->contextHdl)		{ clReleaseContext(pOCL->contextHdl);		pOCL->contextHdl= NULL;	}

	if (bExit)
		exit (iExitCode);
}

void parseArgs(OCLResources* pOCL, int argc, char** argv, unsigned int* test_iterations, char* pDeviceStr, char* pVendorStr, unsigned int* widthReSz, unsigned int* heightReSz, bool* pbShowCL)
{	
	char*			pDeviceWStr = NULL;
	char*			pVendorWStr = NULL;
	const char sUsageString[512] = "Usage: Quicksort [num test iterations] [cpu|gpu] [intel|amd|nvidia] [SurfWidth(^2 only)] [SurfHeight(^2 only)] [show_CL | no_show_CL]";
	
	if (argc != 7)
	{
		Cleanup (pOCL, -1, true, (char *) sUsageString);
	}
	else
	{
		*test_iterations	= atoi (argv[1]);
		pDeviceWStr			= argv[2];			// "cpu" or "gpu"	
		pVendorWStr			= argv[3];			// "intel" or "amd" or "nvidia"
		*widthReSz	= atoi (argv[4]);
		*heightReSz	= atoi (argv[5]);
		if (argv[6][0]=='s')
			*pbShowCL = true;
		else
			*pbShowCL = false;
	}
	sprintf (pDeviceStr, "%s", pDeviceWStr);
	sprintf (pVendorStr, "%s", pVendorWStr);
}

void InstantiateOpenCLKernels(OCLResources *pOCL)
{	
	int t = -1;   // initialize test count (t)

	pOCL->numTests =	++t;

	// Instantiate kernels:

	gqsort_kernel = clCreateKernel(pOCL->programHdl, "gqsort_kernel", &ciErrNum);
	CheckCLError (ciErrNum, "Kernel creation failed.", "Kernel created.");
	lqsort_kernel = clCreateKernel(pOCL->programHdl, "lqsort_kernel", &ciErrNum);
	CheckCLError (ciErrNum, "Kernel creation failed.", "Kernel created.");
}


//#define GET_DETAILED_PERFORMANCE 1
#define RUN_CPU_SORTS
#define HOST 1
#include "Quicksort.h"

uint randomval();

template <class T>
T* partition(T* left, T* right, T pivot) {
    // move pivot to the end
    T temp = *right;
    *right = pivot;
    *left = temp;

    T* store = left;

    for(T* p = left; p != right; p++) {
        if (*p < pivot) {
            temp = *store;
            *store = *p;
            *p = temp;
            store++;
        }
    }

    temp = *store;
    *store = pivot;
    *right = temp;

    return store;
}

template <class T>
void quicksort(T* data, int left, int right)
{
    T* store = partition(data + left, data + right, data[left]);
    int nright = store-data;
    int nleft = nright+1;

    if (left < nright) {
        quicksort(data, left, nright);
	}

    if (nleft < right) {
		quicksort(data, nleft, right); 
	}
}

template <class T>
void gqsort(OCLResources *pOCL, std::vector<block_record>& blocks, std::vector<parent_record>& parents, std::vector<work_record>& news, bool reset) {
	news.resize(blocks.size()*2);

	size_t		dimNDR[2] = { 0, 0};
	size_t		dimWG[2] = { 0, 0 };

	// Create buffer objects for memory.
	cl_mem blocksb = clCreateBuffer(pOCL->contextHdl, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(block_record)*blocks.size(), &blocks[0], &ciErrNum);
	CheckCLError (ciErrNum, "clCreateBuffer failed.", "clCreateBuffer.");
	cl_mem parentsb = clCreateBuffer(pOCL->contextHdl, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(parent_record)*parents.size(), &parents[0], &ciErrNum);
	CheckCLError (ciErrNum, "clCreateBuffer failed.", "clCreateBuffer.");
	cl_mem newsb = clCreateBuffer(pOCL->contextHdl, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(work_record)*news.size(), &news[0], &ciErrNum);
	CheckCLError (ciErrNum, "clCreateBuffer failed.", "clCreateBuffer.");
	

	ciErrNum |= clSetKernelArg(gqsort_kernel, 2, sizeof(cl_mem), (void*) &blocksb);
	CheckCLError(ciErrNum, "clSetKernelArg failed.", "clSetKernelArg");
	ciErrNum |= clSetKernelArg(gqsort_kernel, 3, sizeof(cl_mem), (void*) &parentsb);
	CheckCLError(ciErrNum, "clSetKernelArg failed.", "clSetKernelArg");
	ciErrNum |= clSetKernelArg(gqsort_kernel, 4, sizeof(cl_mem), (void*) &newsb);
	CheckCLError(ciErrNum, "clSetKernelArg failed.", "clSetKernelArg");

	//std::cout << "blocks size is " << blocks.size() << std::endl;
#ifdef GET_DETAILED_PERFORMANCE
	static double absoluteTotal = 0.0;
	static uint count = 0;

	if (reset) {
		absoluteTotal = 0.0;
		count = 0;
	}

	LARGE_INTEGER beginClock, endClock, clockFreq;
	QueryPerformanceFrequency (&clockFreq);
	QueryPerformanceCounter (&beginClock);
#endif
	// Lets do phase 1 pass
	dimNDR[0] = GQSORT_LOCAL_WORKGROUP_SIZE * blocks.size();
	dimWG[0] = GQSORT_LOCAL_WORKGROUP_SIZE;

	ciErrNum = clEnqueueNDRangeKernel (pOCL->cmdQHdl, gqsort_kernel, 1, NULL, dimNDR, dimWG, 0, NULL, 0);
	CheckCLError(ciErrNum, "clEnqueueNDRangeKernel failed.", "clEnqueueNDRangeKernel");
	ciErrNum = clEnqueueReadBuffer(pOCL->cmdQHdl, newsb, CL_TRUE, 0, sizeof(work_record)*news.size(), &news[0], 0, NULL, NULL);
	CheckCLError(ciErrNum, "clEnqueueReadBuffer failed.", "clEnqueueReadBuffer");

#ifdef GET_DETAILED_PERFORMANCE
	QueryPerformanceCounter (&endClock);
	double totalTime = double(endClock.QuadPart - beginClock.QuadPart) / clockFreq.QuadPart;
	absoluteTotal += totalTime;
	std::cout << ++count << ": gqsort time " << absoluteTotal * 1000 << " ms" << std::endl;
#endif
	clReleaseMemObject(blocksb);
	clReleaseMemObject(parentsb);
	clReleaseMemObject(newsb);
}

template <class T>
void lqsort(OCLResources *pOCL, std::vector<work_record>& done, cl_mem tempb, T* d, size_t size) {
	size_t		dimNDR[2] = { 0, 0};
	size_t		dimWG[2] = { 0, 0 };

	//std::cout << "done size is " << done.size() << std::endl; 
	cl_mem doneb = clCreateBuffer(pOCL->contextHdl, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(work_record)*done.size(), &done[0], &ciErrNum);
	CheckCLError (ciErrNum, "clCreateBuffer failed.", "clCreateBuffer.");
	
	ciErrNum |= clSetKernelArg(lqsort_kernel, 2, sizeof(cl_mem), (void*) &doneb);
	CheckCLError(ciErrNum, "clSetKernelArg failed.", "clSetKernelArg");

#ifdef GET_DETAILED_PERFORMANCE
	LARGE_INTEGER beginClock, endClock, clockFreq;
	QueryPerformanceFrequency (&clockFreq);
	QueryPerformanceCounter (&beginClock);
#endif
	// Lets do phase 2 pass
	dimNDR[0] = LQSORT_LOCAL_WORKGROUP_SIZE * (done.size());
	dimWG[0] = LQSORT_LOCAL_WORKGROUP_SIZE;
	ciErrNum = clEnqueueNDRangeKernel (pOCL->cmdQHdl, lqsort_kernel, 1, NULL, dimNDR, dimWG, 0, NULL, 0);
	CheckCLError(ciErrNum, "clEnqueueNDRangeKernel failed.", "clEnqueueNDRangeKernel");

	T* foo = (T*)clEnqueueMapBuffer(pOCL->cmdQHdl, tempb, CL_TRUE, CL_MAP_READ, 0, sizeof(T)*size, 0, 0, 0, &ciErrNum); 

	ciErrNum = clEnqueueUnmapMemObject(pOCL->cmdQHdl, tempb, foo, 0, 0, 0);
	CheckCLError(ciErrNum, "clEnqueueUnmapMemObject failed.", "clEnqueueUnmapMemObject");

#ifdef GET_DETAILED_PERFORMANCE
	QueryPerformanceCounter (&endClock);
	double totalTime = double(endClock.QuadPart - beginClock.QuadPart) / clockFreq.QuadPart;
	std::cout << "lqsort time " << totalTime * 1000 << " ms" << std::endl;
#endif
	clReleaseMemObject(doneb);
}

size_t optp(size_t s, double k, size_t m) {
	return (size_t)pow(2, floor(log(s*k + m)/log(2.0) + 0.5));
}

template <class T>
void GPUQSort(OCLResources *pOCL, size_t size, T* d, T* dn)  {
	// allocate buffers
	cl_mem db = clCreateBuffer(pOCL->contextHdl, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, ((sizeof(T)*size)/64 + 1)*64, d, &ciErrNum);
	CheckCLError (ciErrNum, "clCreateBuffer failed.", "clCreateBuffer.");
	cl_mem dnb = clCreateBuffer(pOCL->contextHdl, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, ((sizeof(T)*size)/64 + 1)*64, dn, &ciErrNum);
	CheckCLError (ciErrNum, "clCreateBuffer failed.", "clCreateBuffer.");

	ciErrNum |= clSetKernelArg(gqsort_kernel, 0, sizeof(cl_mem), (void*) &db);
	CheckCLError(ciErrNum, "clSetKernelArg failed.", "clSetKernelArg");
	ciErrNum |= clSetKernelArg(gqsort_kernel,	1, sizeof(cl_mem), (void*) &dnb);
	CheckCLError(ciErrNum, "clSetKernelArg failed.", "clSetKernelArg");

	ciErrNum |= clSetKernelArg(lqsort_kernel, 0, sizeof(cl_mem), (void*) &db);
	CheckCLError(ciErrNum, "clSetKernelArg failed.", "clSetKernelArg");
	ciErrNum |= clSetKernelArg(lqsort_kernel,	1, sizeof(cl_mem), (void*) &dnb);
	CheckCLError(ciErrNum, "clSetKernelArg failed.", "clSetKernelArg");

	const size_t MAXSEQ = optp(size, 0.00009516, 203);
	const size_t MAX_SIZE = 12*max(MAXSEQ, QUICKSORT_BLOCK_SIZE);
	//std::cout << "MAXSEQ = " << MAXSEQ << std::endl;
	uint startpivot = median(d[0], d[size/2], d[size-1]);
	std::vector<work_record> work, done, news;
	work.reserve(MAX_SIZE);
	done.reserve(MAX_SIZE);
	news.reserve(MAX_SIZE);
	std::vector<parent_record> parent_records;
	parent_records.reserve(MAX_SIZE);
	std::vector<block_record> blocks;
	blocks.reserve(MAX_SIZE);
	
	work.push_back(work_record(0, size, startpivot, 1));

	bool reset = true;

	while(!work.empty() /*&& work.size() + done.size() < MAXSEQ*/) {
		size_t blocksize = 0;
		
		for(auto it = work.begin(); it != work.end(); ++it) {
			blocksize += max((it->end - it->start)/MAXSEQ, 1);
		}
		for(auto it = work.begin(); it != work.end(); ++it) {
			uint start = it->start;
			uint end   = it->end;
			uint pivot = it->pivot;
			uint direction = it->direction;
			uint blockcount = (end - start + blocksize - 1)/blocksize;
			parent_record prnt(start, end, start, end, blockcount-1);
			parent_records.push_back(prnt);

			for(uint i = 0; i < blockcount - 1; i++) {
				uint bstart = start + blocksize*i;
				block_record br(bstart, bstart+blocksize, pivot, direction, parent_records.size()-1);
				blocks.push_back(br);
			}
			block_record br(start + blocksize*(blockcount - 1), end, pivot, direction, parent_records.size()-1);
			blocks.push_back(br);
		}

		gqsort<T>(pOCL, blocks, parent_records, news, reset);
		reset = false;
		//std::cout << " blocks = " << blocks.size() << " parent records = " << parent_records.size() << " news = " << news.size() << std::endl;
		work.clear();
		parent_records.clear();
		blocks.clear();
		for(auto it = news.begin(); it != news.end(); ++it) {
			if (it->direction != EMPTY_RECORD) {
				if (it->end - it->start <= QUICKSORT_BLOCK_SIZE /*size/MAXSEQ*/) {
					if (it->end - it->start > 0)
						done.push_back(*it);
				} else {
					work.push_back(*it);
				}
			}
		}
		news.clear();
	}
	for(auto it = work.begin(); it != work.end(); ++it) {
		if (it->end - it->start > 0)
			done.push_back(*it);
	}

	lqsort<T>(pOCL, done, db, d, size);

	// release buffers: we are done
	clReleaseMemObject(db);
	clReleaseMemObject(dnb);
}

int main(int argc, char** argv)
{
	OCLResources	myOCL;
	unsigned int	test_iterations;
	char			pDeviceStr[256];
	char			pVendorStr[256];
	const char*		pSourceFileStr	= "QuicksortKernels.cl";
	bool			bShowCL = false;

	// Image data:
	uint			heightReSz, widthReSz;

	double totalTime, quickSortTime, stdSortTime;

	LARGE_INTEGER beginClock, endClock, clockFreq;
	QueryPerformanceFrequency (&clockFreq);

	parseArgs (&myOCL, argc, argv, &test_iterations, pDeviceStr, pVendorStr, &widthReSz, &heightReSz, &bShowCL);

	printf("\n\n\n--------------------------------------------------------------------\n");
	
	uint arraySize = widthReSz*heightReSz;
	printf("Allocating array size of %d\n", arraySize);
	uint* pArray = (uint*)_aligned_malloc (((arraySize*sizeof(uint))/64 + 1)*64, 4096);
	uint* pArrayCopy = (uint*)_aligned_malloc (((arraySize*sizeof(uint))/64 + 1)*64, 4096);

	std::generate(pArray, pArray + arraySize, []() { static uint j = 0; return j++; });
	std::random_shuffle(pArray, pArray + arraySize);
#ifdef RUN_CPU_SORTS
	std::cout << "Sorting the regular way..." << std::endl;
	std::copy(pArray, pArray + arraySize, pArrayCopy);

	QueryPerformanceCounter (&beginClock);
	std::sort(pArrayCopy, pArrayCopy + arraySize);
	QueryPerformanceCounter (&endClock);
	totalTime = double(endClock.QuadPart - beginClock.QuadPart) / clockFreq.QuadPart;	
	std::cout << "Time to sort: " << totalTime * 1000 << " ms" << std::endl;
	stdSortTime = totalTime;

	std::cout << "Sorting with quicksort on the cpu: " << std::endl;
	std::copy(pArray, pArray + arraySize, pArrayCopy);

	QueryPerformanceCounter (&beginClock);
	quicksort(pArrayCopy, 0, arraySize-1);
	QueryPerformanceCounter (&endClock);
	totalTime = double(endClock.QuadPart - beginClock.QuadPart) / clockFreq.QuadPart;
	std::cout << "Time to sort: " << totalTime * 1000 << " ms" << std::endl;
	quickSortTime = totalTime;
#ifdef TRUST_BUT_VERIFY
	{
		std::vector<uint> verify(arraySize);
		std::copy(pArray, pArray + arraySize, verify.begin());

		std::cout << "verifying: ";
		std::sort(verify.begin(), verify.end());
		bool correct = std::equal(verify.begin(), verify.end(), pArrayCopy);
		if (!correct) {
			for(size_t i = 0; i < arraySize; i++) {
				if (verify[i] != pArrayCopy[i]) {
					std:: cout << "discrepancy at " << i << " " << pArrayCopy[i] << std::endl;
				}
			}
		}
		std::cout << std::boolalpha << correct << std::endl;
	}
#endif
#endif // RUN_CPU_SORTS

	// Initialize OpenCL:
	InitializeOpenCL (pDeviceStr, pVendorStr, &myOCL.deviceID, &myOCL.contextHdl, &myOCL.cmdQHdl);
	if (bShowCL)
		QueryPrintOpenCLDeviceInfo (myOCL.deviceID, myOCL.contextHdl);	
	QueryPerformanceCounter (&beginClock);
	CompileOpenCLProgram (myOCL.deviceID, myOCL.contextHdl, pSourceFileStr, &myOCL.programHdl);
	QueryPerformanceCounter (&endClock);
	totalTime = double(endClock.QuadPart - beginClock.QuadPart) / clockFreq.QuadPart;
	std::cout << "Time to build OpenCL Program: " << totalTime * 1000 << " ms" << std::endl;
	InstantiateOpenCLKernels (&myOCL);

	std::cout << "Sorting with GPUQSort on the " << pDeviceStr << ": " << std::endl;

	std::vector<uint> original(arraySize);
	std::copy(pArray, pArray + arraySize, original.begin());

	uint NUM_ITERATIONS = test_iterations;
	std::vector<double> times;
	times.resize(NUM_ITERATIONS);
	double AverageTime = 0.0;
	for(uint k = 0; k < NUM_ITERATIONS; k++) {
		std::copy(original.begin(), original.end(), pArray);
		//std::copy(pArray, pArray + arraySize, pArrayCopy);
		std::vector<uint> seqs;
		std::vector<uint> verify(arraySize);
		std::copy(pArray, pArray + arraySize, verify.begin());

		QueryPerformanceCounter (&beginClock);
		GPUQSort(&myOCL, arraySize, pArray, pArrayCopy);
		QueryPerformanceCounter (&endClock);
		totalTime = double(endClock.QuadPart - beginClock.QuadPart) / clockFreq.QuadPart;
		std::cout << "Time to sort: " << totalTime * 1000 << " ms" << std::endl;
		times[k] = totalTime;
		AverageTime += totalTime;
#ifdef TRUST_BUT_VERIFY
		std::cout << "verifying: ";
		std::sort(verify.begin(), verify.end());
		bool correct = std::equal(verify.begin(), verify.end(), pArray);
		if (!correct) {
			for(size_t i = 0; i < arraySize; i++) {
				if (verify[i] != pArray[i]) {
					std:: cout << "discrepancy at " << i << " " << pArray[i] << std::endl;
				}
			}
		}
		std::cout << std::boolalpha << correct << std::endl;
#endif
	}
	AverageTime = AverageTime/NUM_ITERATIONS;
	std::cout << "Average Time: " << AverageTime * 1000 << " ms" << std::endl;
	double stdDev = 0.0, minTime = 1000000.0, maxTime = 0.0;
	for(uint k = 0; k < NUM_ITERATIONS; k++) 
	{
		stdDev += (AverageTime - times[k])*(AverageTime - times[k]);
		minTime = min(minTime, times[k]);
		maxTime = max(maxTime, times[k]);
	}

	stdDev = sqrt(stdDev/(NUM_ITERATIONS - 1));
	std::cout << "Standard Deviation: " << stdDev * 1000 << std::endl;
	std::cout << "%error (3*stdDev)/Average: " << 3*stdDev / AverageTime * 100 << "%" << std::endl;
	std::cout << "min time: " << minTime * 1000 << " ms" << std::endl;
	std::cout << "max time: " << maxTime * 1000 << " ms" << std::endl;
#ifdef RUN_CPU_SORTS
	std::cout << "Average speedup over CPU quicksort: " << quickSortTime/AverageTime << std::endl;
	std::cout << "Average speedup over CPU std::sort: " << stdSortTime/AverageTime << std::endl;
#endif // RUN_CPU_SORTS

	printf("-------done--------------------------------------------------------\n");
	getchar();
	Sleep(2000);

	_aligned_free(pArray);
	_aligned_free(pArrayCopy);

	return 0;
}
