#include <microtime.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <vector>
#include <time.h>
#include <mpi.h>
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

// -------------------------------------------------------------------------------------------//

int main(int argc, char** argv) {

	MPI_Init(&argc, &argv);

	int comm_size;
	int proc_num;

	MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &proc_num);

	if (argc != 4) {
		if (proc_num == 0)
			fprintf(stderr, "USAGE: %s <numNums> <randomNumBound> <samplesPerBucket>\n", argv[0]);
		MPI_Finalize();
		exit(1);
	}

//	double t, time1, time2;

//	char strBuff0[10000] = "      ";

//	char strBuff[10000] =  "      ";

	int inputArraySize = atoi(argv[1]);
	int randomNumBound = atoi(argv[2]);
	int samplesPerBucket = atoi(argv[3]);

//	int numsTruthSort[inputArraySize];
	int inputSS[inputArraySize];

	int numNums = sizeof(inputSS)/sizeof(inputSS[0]);

//-------------------------------------------//

	if (proc_num == 0) {
		generateRandNums(inputSS, inputArraySize, randomNumBound);
/*COPY ARRAY*/
//		memcpy(numsTruthSort, inputSS, numNums * sizeof(int));

		//     sprintf(strBuff0 + strlen(strBuff0),"\nRandomly Generated Numbers:\n");
		//     printArray(numsTruthSort, numNums, strBuff0);
	}

//-------------------------------------------//

int prepState = -2;
int splitters[comm_size-1];
std::vector<int> bucket;

	if (proc_num == 0) {

/*SERIAL BUBBLE SORT*/
//		bubblesort(numsTruthSort,0,numNums);

/*SERIAL QUICK SORT*/
//		quicksort(numsTruthSort,0,numNums - 1);

		//     sprintf(strBuff0 + strlen(strBuff0), "\nSorted Numbers (using serial bubble sort):\n");
		//     printArray(numsTruthSort, numNums, strBuff0);

		prepState = samplesort_prep(inputSS, numNums, comm_size, samplesPerBucket);

			if (prepState == 1)
			samplesort_get_splitters(inputSS, numNums, comm_size, samplesPerBucket, splitters);

	}

	MPI_Bcast(&prepState,1,MPI_INT,0,MPI_COMM_WORLD);

	if (prepState == 1) {

		MPI_Bcast(inputSS,numNums,MPI_INT,0,MPI_COMM_WORLD);
		MPI_Bcast(splitters,comm_size-1,MPI_INT,0,MPI_COMM_WORLD);

		samplesort_get_sorted_buckets(proc_num, inputSS, numNums, comm_size, samplesPerBucket, splitters, bucket);

		//     sprintf(strBuff + strlen(strBuff), "\nSorted Bucket (proc # %i):\n", proc_num);
		//     printArray(&bucket[0],bucket.size(), strBuff);

		int bucketSize = bucket.size();
		int bucket_sizes[comm_size];

		//     sprintf(strBuff + strlen(strBuff), "\nBucket size:%i\n", bucketSize);

		MPI_Gather(&bucketSize,1,MPI_INT,bucket_sizes,1,MPI_INT,0,MPI_COMM_WORLD);

		int displc[comm_size];

		if (proc_num == 0) {
			displc[0] = 0;
			for (int i = 1; i < comm_size; ++i) {
				displc[i] = displc[i-1] + bucket_sizes[i-1];
			}
		}

		MPI_Gatherv(bucket.data(),bucketSize,MPI_INT,inputSS,bucket_sizes,displc,MPI_INT,0,MPI_COMM_WORLD);
	}

//-------------------------------------------//
/*
	int displc[comm_size];
	int str_buffer_sizes[comm_size];
	int len_str_buffer = strlen(strBuff);

	MPI_Gather(&len_str_buffer,1,MPI_INT,str_buffer_sizes,1,MPI_INT,0,MPI_COMM_WORLD);

	int outputStringSize;
	char *outputString;

	if (proc_num == 0) {
		displc[0] = 0;
		outputStringSize = str_buffer_sizes[0];
		for (int i = 1; i < comm_size; ++i) {
			outputStringSize += str_buffer_sizes[i];
			displc[i] = displc[i-1] + str_buffer_sizes[i-1];
		}
		outputString = (char *)malloc((outputStringSize+1) * sizeof(char));
	}

	MPI_Gatherv(&strBuff,strlen(strBuff),MPI_CHAR,outputString,str_buffer_sizes,displc,MPI_CHAR,0,MPI_COMM_WORLD);
*/
	if (proc_num == 0) {

		//     sprintf(strBuff0 + strlen(strBuff0), "\nSorted Numbers (using parallel sample sort sort):\n");
		//     printArray(inputSS, numNums, strBuff0);

/*COMPARE*/
//		int sum_abs_diff = compareArray(numsTruthSort,numNums, inputSS, numNums);

/*SPRINTF*/
//		sprintf(strBuff0 + strlen(strBuff0), "sum of abs value of diff = %d \n\n",sum_abs_diff);

		printf("\n------------------------------------------------------------------\n");
//		printf("%s\n",strBuff0);
/* CALC EVEN-ODD DIFF */
		printf("Problem size: %i, # of processes: %i\n", numNums, comm_size);
		printf("Even Odd Diff Calc: %i", calcEvenOddDiff(inputSS, numNums));
		printf("\n------------------------------------------------------------------\n");
//		printf("%s\n",outputString);
	}

//	free(outputString);
	MPI_Finalize();
	return 0;
}
