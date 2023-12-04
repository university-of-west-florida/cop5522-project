#include <microtime.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <vector>
#include <time.h>
#include <omp.h>
#include <string.h>

// -------------------------------------------------------------------------------------------//

void printArray(int arrayIn[], int arraySize, char strBuff[]) {
    sprintf(strBuff + strlen(strBuff), "Array Values:\n\t");
    for(int i = 0; i < arraySize; ++i)
        sprintf (strBuff + strlen(strBuff), "%d ",arrayIn[i]);
    sprintf(strBuff + strlen(strBuff), "\n\n");
}


// -------------------------------------------------------------------------------------------//

void printPartialArray(int arrayIn[], int start, int end, char strBuff[]) {
    sprintf(strBuff + strlen(strBuff), "Partial Array Values (i=%i - %i):\n\t",start,end);
    for(int i = start; i < end; ++i)
        sprintf (strBuff + strlen(strBuff), "%d ",arrayIn[i]);
    sprintf(strBuff + strlen(strBuff), "\n\n");
}


// -------------------------------------------------------------------------------------------//

int compareArray(int array1[], int arraySize1, int array2[], int arraySize2) {
	if (arraySize1 != arraySize2)
		return -30303030;
	int sum_abs_diff = 0;
	for(int i = 0; i < arraySize1; ++i)
		sum_abs_diff += abs(array1[i]-array2[i]);
	return sum_abs_diff;
}

// -------------------------------------------------------------------------------------------//

int calcEvenOddDiff(int arrayIn[], int arraySize) {
	int calc = 0;
	for (int i = 0; i < arraySize ; ++i) {
		if (i % 2 == 0)
			calc += arrayIn[i];
		else
			calc -= arrayIn[i];
	}
	return calc;
}
// -------------------------------------------------------------------------------------------//

/*
void generateRandNums(int numbersToKeep[], int numbersToSort[], int numNums, int randomNumBound) {
        int i;
   	srand(1);
	for(i = 0; i < numNums; ++i) {
		numbersToKeep[i] = rand() % randomNumBound;
		numbersToSort[i] = numbersToKeep[i];
	}
}
*/

void generateRandNums(int arrayRandNums[], int numNums, int randomNumBound) {
        int i;
   	srand(1);
	for(i = 0; i < numNums; ++i) {
		arrayRandNums[i] = rand() % randomNumBound;
	}
}

// -------------------------------------------------------------------------------------------//

void bubblesort(int inputNums[], int start, int end) {
    int i,j,temp;
    for(i = start; i < end - 1; ++i) {
        for(j = start; j < start + end - i - 1; ++j) {
            if(inputNums[j] > inputNums[j + 1]) {
                temp = inputNums[j + 1];
                inputNums[j + 1] = inputNums[j];
                inputNums[j] = temp;
            }
        }
    }
}

// -------------------------------------------------------------------------------------------//

// from geeksforgeeks.org/quick-sort/

int partition(int arr[],int low,int high)
{
  int pivot=arr[high];
  int i=(low-1);

  for(int j=low;j<=high;j++)
  {
    if(arr[j]<pivot)
    {
      i++;
      std::swap(arr[i],arr[j]);
    }
  }
  std::swap(arr[i+1],arr[high]);
  return (i+1);
}

void quicksort(int arr[],int low,int high)
{
  if(low<high)
  {

    int p_i=partition(arr,low,high);

    quicksort(arr,low,p_i-1);
    quicksort(arr,p_i+1,high);
  }
}

// -------------------------------------------------------------------------------------------//

int samplesort_prep(int arrayIn[], int arraySize, int numOfBuckets, int samplesPerBucket) {
	if (numOfBuckets == 1 || numOfBuckets >= arraySize)  {
		quicksort(arrayIn, 0, arraySize-1);
		return -1;
	}
	return 1;
}

// -------------------------------------------------------------------------------------------//

void samplesort_get_splitters(int arrayIn[], int arraySize, int numOfBuckets, int samplesPerBucket, int splittersOut[]) {

	int i;

	if(numOfBuckets * samplesPerBucket > arraySize)
		samplesPerBucket = arraySize / numOfBuckets;

	int sortedSamples[(numOfBuckets - 1) * samplesPerBucket];

	for(i = 0; i < (numOfBuckets - 1) * samplesPerBucket; ++i)
		sortedSamples[i] = arrayIn[i];

	quicksort(sortedSamples, 0, (numOfBuckets - 1) * samplesPerBucket - 1);

	for(i = 0; i < numOfBuckets - 1; ++i)
		splittersOut[i] = sortedSamples[i * samplesPerBucket];
}

// -------------------------------------------------------------------------------------------//
/*
void samplesort_get_sorted_buckets(int proc_num, int arrayIn[], int arraySize, int numOfBuckets, int samplesPerBucket, int splitters[], std::vector<int>& bucket) {

	int i;

	if (proc_num == 0) {
		for (i = 0; i < arraySize; ++i) {
			if (arrayIn[i] <= splitters[proc_num]) {
				bucket.push_back(arrayIn[i]);
			}
		}
	}
	else if(proc_num < numOfBuckets - 1) {
		for (i = 0; i <= arraySize; ++i) {
			if ((arrayIn[i] <= splitters[proc_num]) && (arrayIn[i] > splitters[proc_num-1]))
				bucket.push_back(arrayIn[i]);
		}
	}
	else {
		for (i = 0; i < arraySize; ++i) {
			if (arrayIn[i] > splitters[proc_num - 1])
				bucket.push_back(arrayIn[i]);
		}
	}

	quicksort(&bucket[0],0,bucket.size()-1);
}
*/
// -------------------------------------------------------------------------------------------//

int main(int argc, char** argv) {

	if (argc != 5) {
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

	int* arrayIn = (int*) malloc (arraySize * sizeof(int));

//-------------------------------------------//

//	generateRandNums(inputSS, inputArraySize, randomNumBound);
	generateRandNums(arrayIn, arraySize, randomNumBound);

/*COPY ARRAY*/
//		memcpy(numsTruthSort, inputSS, numNums * sizeof(int));

		//     sprintf(strBuff0 + strlen(strBuff0),"\nRandomly Generated Numbers:\n");
		//     printArray(numsTruthSort, numNums, strBuff0);

//-------------------------------------------//

int prepState = -2;
int splitters[comm_size-1];
int bucketSizes[numOfBuckets] = {0};
int dataLocations[numOfBuckets] = {0};

/*SERIAL BUBBLE SORT*/
//		bubblesort(numsTruthSort,0,numNums);

/*SERIAL QUICK SORT*/
//		quicksort(numsTruthSort,0,numNums - 1);

		//     sprintf(strBuff0 + strlen(strBuff0), "\nSorted Numbers (using serial bubble sort):\n");
		//     printArray(numsTruthSort, numNums, strBuff0);

	prepState = samplesort_prep(arrayIn, arraySize, comm_size, samplesPerBucket);

	if (prepState == 1) {

		samplesort_get_splitters(arrayIn, arraySize, comm_size, samplesPerBucket, splitters);

		int i;
		std::vector<int> bucket;
		int proc_num;

		#pragma omp parallel \
			shared(arrayIn, arraySize, splitters, numOfBuckets, bucketSizes, dataLocations) private(i, bucket, proc_num)
		{
			proc_num = omp_get_thread_num();

			if (proc_num == 0) {
				for (i = 0; i < arraySize; ++i) {
					if (arrayIn[i] <= splitters[proc_num]) {
						bucket.push_back(arrayIn[i]);
					}
				}
			}
			else if(proc_num < numOfBuckets - 1) {
				for (i = 0; i <= arraySize; ++i) {
					if ((arrayIn[i] <= splitters[proc_num]) && (arrayIn[i] > splitters[proc_num-1]))
						bucket.push_back(arrayIn[i]);
				}
			}
			else {
				for (i = 0; i < arraySize; ++i) {
					if (arrayIn[i] > splitters[proc_num - 1])
						bucket.push_back(arrayIn[i]);
				}
			}

		       quicksort(&bucket[0],0,bucket.size()-1);

			#pragma omp barrier

			#pragma omp critical
			{
				bucketSizes[proc_num] = bucket.size();
			}

			#pragma omp barrier

			if (proc_num == 0) {
				for (i = 1; i < numOfBuckets; ++i) {
					dataLocations[i] = dataLocations[i-1] + bucketSizes[i-1];
				}
			}

			#pragma omp barrier

			#pragma omp critical
			{
			    memcpy(&arrayIn[dataLocations[proc_num]],bucket.data(),bucket.size() * sizeof(int));
			}
		}

	}

//-------------------------------------------//

	printf("\n------------------------------------------------------------------\n");

        printf("Problem size: %i, # of processes: %i\t\t", arraySize, comm_size);
        printf("(Random number bound: %i, samplesPerBucket: %i)\n", randomNumBound, samplesPerBucket);
        printf("Even Odd Diff Calc: %i\t\t\t\t", calcEvenOddDiff(arrayIn, arraySize));

	free(arrayIn);

        time2 = microtime();
        t = time2 - time1;

        printf("(Process 0 time: %0.3f sec)\n", t * 1e-6);

	return 0;
}
