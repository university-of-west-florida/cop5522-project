#pragma once
#include "microtime.h"

#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <unistd.h>
#include <vector>
#include <time.h>
#include <pthread.h>
#include <typeinfo>

#define DEBUG(...) currentDebugNames = #__VA_ARGS__; printDebug(__VA_ARGS__)
std::string currentDebugNames;
#define ENDL std::cout << std::endl

#define INDICES_PER_THREAD 10

// Holding this for potentially threading out the entire sort algorithm
struct sampleSortArgs
{
	std::vector<int>* inputNums;
	std::vector<pthread_t*>* threads;
	int start;
	int end;
	int numOfBuckets;
	int samplesPerBucket;
	pthread_mutex_t mtx;
};

// Will rename depending on if the above is kept
struct testArgs
{
	std::vector<int>* inputNums;
	std::vector<int>* splitters;
	std::vector<std::vector<int>>* buckets;
	int i, j, k, start, numNums, numOfBuckets;
	pthread_mutex_t mtx;
};
