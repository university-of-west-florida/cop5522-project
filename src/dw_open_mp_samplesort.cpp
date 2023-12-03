#include <microtime.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <vector>
#include <time.h>
#include <omp.h>
#include <string.h>

// -------------------------------------------------------------------------------------------//
// Print the elements of an array along with a given string buffer.
void printArray(int arrayIn[], int arraySize, char strBuff[])
{
	sprintf(strBuff + strlen(strBuff), "Array Values:\n\t");
	for (int i = 0; i < arraySize; ++i)
		sprintf(strBuff + strlen(strBuff), "%d ", arrayIn[i]);
	sprintf(strBuff + strlen(strBuff), "\n\n");
}

// -------------------------------------------------------------------------------------------//
// Print specified range of elements from an array along with a given string buffer.
void printPartialArray(int arrayIn[], int start, int end, char strBuff[])
{
	sprintf(strBuff + strlen(strBuff), "Partial Array Values (i=%i - %i):\n\t", start, end);
	for (int i = start; i < end; ++i)
		sprintf(strBuff + strlen(strBuff), "%d ", arrayIn[i]);
	sprintf(strBuff + strlen(strBuff), "\n\n");
}

// -------------------------------------------------------------------------------------------//
// Comparing two arrays and return the sum of absolute differences between corresponding elements.
int compareArray(int array1[], int arraySize1, int array2[], int arraySize2)
{
	if (arraySize1 != arraySize2)
		return -30303030;
	int sum_abs_diff = 0;
	for (int i = 0; i < arraySize1; ++i)
		sum_abs_diff += abs(array1[i] - array2[i]);
	return sum_abs_diff;
}

// -------------------------------------------------------------------------------------------//
// Calculate the difference between the sum of even-indexed and odd-indexed elements in an array.
int calcEvenOddDiff(int arrayIn[], int arraySize)
{
	int calc = 0;
	for (int i = 0; i < arraySize; ++i)
	{
		if (i % 2 == 0)
			calc += arrayIn[i];
		else
			calc -= arrayIn[i];
	}
	return calc;
}

// -------------------------------------------------------------------------------------------//
// Checking if an array is sorted in ascending order.
bool sort_check(int arrayIn[], int arraySize)
{
	for (int i = 1; i < arraySize; ++i)
	{
		if (arrayIn[i - 1] > arrayIn[i])
			return false;
	}
	return true;
}

// -------------------------------------------------------------------------------------------//
// Generate an array of random numbers within a specified range.
void generateRandNums(int arrayRandNums[], int numNums, int randomNumBound)
{
	int i;
	srand(1);
	for (i = 0; i < numNums; ++i)
	{
		arrayRandNums[i] = rand() % randomNumBound;
	}
}

// -------------------------------------------------------------------------------------------//
// Perform the bubble sort algorithm on a specified range of elements in an array.
void bubblesort(int inputNums[], int start, int end)
{
	int i, j, temp;
	for (i = start; i < end - 1; ++i)
	{
		for (j = start; j < start + end - i - 1; ++j)
		{
			if (inputNums[j] > inputNums[j + 1])
			{
				temp = inputNums[j + 1];
				inputNums[j + 1] = inputNums[j];
				inputNums[j] = temp;
			}
		}
	}
}

// -------------------------------------------------------------------------------------------//

// from geeksforgeeks.org/quick-sort/
// Perform the partition step in the quicksort algorithm.
int partition(int arr[], int low, int high)
{
	int pivot = arr[high];
	int i = (low - 1);

	for (int j = low; j <= high; j++)
	{
		if (arr[j] < pivot)
		{
			i++;
			std::swap(arr[i], arr[j]);
		}
	}
	std::swap(arr[i + 1], arr[high]);
	return (i + 1);
}
// Implementing quicksort algorithm on a specified range of elements in an array.
void quicksort(int arr[], int low, int high)
{
	if (low < high)
	{

		int p_i = partition(arr, low, high);

		quicksort(arr, low, p_i - 1);
		quicksort(arr, p_i + 1, high);
	}
}

// -------------------------------------------------------------------------------------------//
// Prepare samplesort algorithm by checking the number of buckets and performing quicksort if necessary.
int samplesort_prep(int arrayIn[], int arraySize, int numOfBuckets, int samplesPerBucket)
{
	if (numOfBuckets == 1 || numOfBuckets >= arraySize)
	{
		quicksort(arrayIn, 0, arraySize - 1);
		return -1;
	}
	return 1;
}

// -------------------------------------------------------------------------------------------//
// Get the splitters for the samplesort algorithm based on the given array, number of buckets, and samples per bucket.
void samplesort_get_splitters(int arrayIn[], int arraySize, int numOfBuckets, int samplesPerBucket, int splittersOut[])
{

	int i;

	if (numOfBuckets * samplesPerBucket > arraySize)
		samplesPerBucket = arraySize / numOfBuckets;

	int sortedSamples[(numOfBuckets - 1) * samplesPerBucket];

	for (i = 0; i < (numOfBuckets - 1) * samplesPerBucket; ++i)
		sortedSamples[i] = arrayIn[i];

	quicksort(sortedSamples, 0, (numOfBuckets - 1) * samplesPerBucket - 1);

	for (i = 0; i < numOfBuckets - 1; ++i)
		splittersOut[i] = sortedSamples[i * samplesPerBucket];
}

// -------------------------------------------------------------------------------------------//
// Distribute array elements into partial buckets based on the splitters.
void make_partial_buckets(int arrayIn[], int arraySize, int numOfBuckets, int splitters[], std::vector<int> B[], int B_sizes[])
{

	int i, k;

	for (i = 0; i < arraySize; ++i)
	{
		if (arrayIn[i] > splitters[numOfBuckets - 2])
		{
			B[numOfBuckets - 1].push_back(arrayIn[i]);
			B_sizes[numOfBuckets - 1] += 1;
		}
		else
		{
			for (k = 0; k < numOfBuckets - 1; ++k)
			{
				if (arrayIn[i] <= splitters[k])
				{
					B[k].push_back(arrayIn[i]);
					B_sizes[k] += 1;
					break;
				}
			}
		}
	}
}

// -------------------------------------------------------------------------------------------//
// Main function to execute the entire program, including array generation, sorting, and timing measurements.
int main(int argc, char **argv)
{

	if (argc != 5)
	{
		fprintf(stderr, "USAGE: %s <numOfProcessors> <numNums> <randomNumBound> <samplesPerBucket>\n", argv[0]);
		exit(1);
	}

	double t, time1, time2;

	time1 = microtime();

	int comm_size = atoi(argv[1]);
	int arraySize = atoi(argv[2]);
	int randomNumBound = atoi(argv[3]);
	int samplesPerBucket = atoi(argv[4]);

	int numOfBuckets = comm_size;

	omp_set_dynamic(0);
	omp_set_num_threads(comm_size);

	int *arrayIn = (int *)malloc(arraySize * sizeof(int));

	//-------------------------------------------//

	generateRandNums(arrayIn, arraySize, randomNumBound);

	//-------------------------------------------//

	int prepState = -2;

	prepState = samplesort_prep(arrayIn, arraySize, comm_size, samplesPerBucket);

	if (prepState == 1)
	{

		int splitters[comm_size - 1];
		samplesort_get_splitters(arrayIn, arraySize, comm_size, samplesPerBucket, splitters);

		int bucketSizes[numOfBuckets] = {0};
		int bucketStartLocations[numOfBuckets] = {0};
		int bucketOffsets[numOfBuckets] = {0};

		int i;
		int proc_num;

		std::vector<int> B[comm_size];
		int B_size[comm_size] = {0};

#pragma omp parallel                                                                                                                     \
shared(arrayIn, arraySize, comm_size, splitters, numOfBuckets, bucketSizes, bucketStartLocations, bucketOffsets) private(i, proc_num, B) \
	firstprivate(B_size)
		{
			proc_num = omp_get_thread_num();

			make_partial_buckets(&arrayIn[proc_num * (arraySize / comm_size)], arraySize / comm_size, numOfBuckets, splitters, B, B_size);

			if (proc_num == 0)
				make_partial_buckets(&arrayIn[arraySize - (arraySize % comm_size)], arraySize % comm_size, numOfBuckets, splitters, B, B_size);

#pragma omp critical
			{
				for (i = 0; i < comm_size; ++i)
					bucketSizes[i] += B_size[i];
			}

#pragma omp barrier

			if (proc_num == 0)
			{
				for (i = 1; i < numOfBuckets; ++i)
				{
					bucketStartLocations[i] = bucketStartLocations[i - 1] + bucketSizes[i - 1];
				}
			}

#pragma omp barrier

#pragma omp critical
			{
				for (i = 0; i < comm_size; ++i)
				{
					memcpy(&arrayIn[bucketStartLocations[i] + bucketOffsets[i]], B[i].data(), B[i].size() * sizeof(int));
					bucketOffsets[i] += B_size[i];
				}
			}

#pragma omp barrier

			quicksort(arrayIn, bucketStartLocations[proc_num], bucketStartLocations[proc_num] + bucketSizes[proc_num] - 1);
		}
	}

	//-------------------------------------------//

	printf("\n------------------------------------------------------------------\n");

	printf("Problem size: %i, # of processes: %i\t\t", arraySize, comm_size);
	printf("(Random number bound: %i, samplesPerBucket: %i)\n", randomNumBound, samplesPerBucket);
	printf("Sort successful: %s\t\t\t\t\t", sort_check(arrayIn, arraySize) ? "TRUE" : "FALSE");
	free(arrayIn);

	time2 = microtime();
	t = time2 - time1;

	printf("(Process 0 time: %0.3f sec)\n", t * 1e-6);

	return 0;
}
