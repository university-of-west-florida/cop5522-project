#include "microtime.h"

#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <unistd.h>
#include <vector>
#include <time.h>
#include <pthread.h>
#include <typeinfo>
#include <errno.h>
#include <string>
#include <string.h>
#include <cstdlib>
#include <algorithm>
#include <random>
#include <unordered_map>


#define DEBUG(...) currentDebugNames = #__VA_ARGS__; printDebug(__VA_ARGS__)
std::string currentDebugNames;
#define ENDL std::cout << std::endl


struct buildBucketsArgs
{
	std::vector<int>* valsToSort;
	std::vector<std::vector<int>>* buckets;
	std::vector<std::vector<int>>* splitters;
	long long offset;
	long long indexRange;
	int numOfBuckets;
	int samplesPerBucket;
	pthread_mutex_t* mtx;
};

struct sortAndCombineArgs
{
	std::vector<int>* bucket;
	std::vector<int>* putHere;
	long long offset;
	long long numToDo;
};

struct generateVectorArgs
{
	std::vector<int>* vec;
	long long offset;
	long long indexRange;
	int randSeed;
	int upperBound;
};
