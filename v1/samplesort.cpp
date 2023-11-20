#include <microtime.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <vector>
#include <time.h>
#include <mpi.h>
#include <string.h>

// -------------------------------------------------------------------------------------------//

void printArray(int inputNums[], int numNums, char strBuff[]) {
    sprintf(strBuff + strlen(strBuff), "Array Values:\n\t");
    int i;
    for(i = 0; i < numNums; ++i)
        sprintf (strBuff + strlen(strBuff), "%d ",inputNums[i]);
    sprintf(strBuff + strlen(strBuff), "\n\n");
}

// -------------------------------------------------------------------------------------------//

int compareArray(int origNumbers[], int sortedNumbers[], int numNums) {
    int i;
    int sum_abs_diff = 0;
    for(i = 0; i < numNums; ++i)
        sum_abs_diff += abs(origNumbers[i]-sortedNumbers[i]);
    return sum_abs_diff;
}

// -------------------------------------------------------------------------------------------//

void generateRandNums(int numbersToKeep[], int numbersToSort[], int numNums, int randomNumBound) {
        int i;
   	srand(1);
	for(i = 0; i < numNums; ++i) {
		numbersToKeep[i] = rand() % randomNumBound;
		numbersToSort[i] = numbersToKeep[i];
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

void getMinMax(int inputNums[], int inputSize, int minMax[]) {
        minMax[0] = inputNums[0];
	minMax[1] = minMax[0];
   	for(int i = 1; i < inputSize; ++i) {
		if (inputNums[i] < minMax[0])
			minMax[0] = inputNums[i];
		if (inputNums[i] > minMax[1])
			minMax[1] = inputNums[i];
	}
}

// -------------------------------------------------------------------------------------------//

void samplesort(int inputNums[], int start, int end, int numOfBuckets, int samplesPerBucket) {
	if (numOfBuckets == 1)  {
		bubblesort(inputNums, start, end);
		return;
	}

	int i, j, k;
        int numNums = end - start;

            if(numOfBuckets * samplesPerBucket > numNums) {
        if(numNums / numOfBuckets >= 1) {
            samplesPerBucket = 1;
        }
                else {
            bubblesort(inputNums, start, end);
            return;
        }
    }

            int sortedSamples[(numOfBuckets - 1) * samplesPerBucket];

	for(i = start; i < start + (numOfBuckets - 1) * samplesPerBucket; ++i)
        	sortedSamples[i - start] = inputNums[i];

        bubblesort(sortedSamples, 0, (numOfBuckets - 1) * samplesPerBucket);

        int splitters[numOfBuckets - 1];
    for(i = 0; i < numOfBuckets - 1; ++i)
        splitters[i] = sortedSamples[i * samplesPerBucket];

                    std::vector<int> buckets[numOfBuckets];
        for(i = 0; i < numOfBuckets; ++i) {
        for(j = start + (i * (numNums / numOfBuckets)); j < start + ((i + 1) * (numNums / numOfBuckets)); ++j) {
                                    for(k = 0; k < numOfBuckets - 1; ++k) {
                if(inputNums[j] <= splitters[k]) {
                    buckets[k].push_back(inputNums[j]);
                    break;
                }
                if(k == numOfBuckets - 2)
                    buckets[numOfBuckets - 1].push_back(inputNums[j]);
            }
        }
    }

                    for(i = start + numNums - (numNums % numOfBuckets); i < end; ++i) {
        for(k = 0; k < numOfBuckets - 1; ++k) {
            if(inputNums[i] <= splitters[k]) {
                buckets[k].push_back(inputNums[i]);
                break;
            }
            if(k == numOfBuckets - 2)
                buckets[numOfBuckets - 1].push_back(inputNums[i]);
        }
    }

            int runningIndex = start;
    for(i = 0; i < numOfBuckets; ++i) {
        int bucketSize = buckets[i].size();
        for(j = 0; j < bucketSize; ++j )
            inputNums[runningIndex + j] = buckets[i].at(j);
                        samplesort(inputNums, runningIndex, runningIndex + bucketSize, numOfBuckets, samplesPerBucket);
        runningIndex += bucketSize;
    }
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
			fprintf(stderr, "USAGE: %s <numNums> <randomNumBound> <num of processes>\n", argv[0]);
		MPI_Finalize();
		exit(1);
	}

	double t, time1, time2;
	char strBuff0[10000];

	char strBuff[10000];

	int inputArraySize = atoi(argv[1]);
	int randomNumBound = atoi(argv[2]);
	int numOfProcesses = atoi(argv[3]);
	int samplesPerBucket = 5;

	int numsTruthSort[inputArraySize];
	int inputSS[inputArraySize];

	int numNums = sizeof(inputSS)/sizeof(inputSS[0]);

	int elements_per_proc[comm_size];
	int displc[comm_size];
 	int rcv_buff_size = numNums/comm_size;

	if (proc_num == comm_size - 1) {
		rcv_buff_size += numNums % comm_size;
	}

	if (proc_num == 0) {
		generateRandNums(numsTruthSort, inputSS, inputArraySize, randomNumBound);
		for (int i = 0; i < comm_size; ++i) {
			elements_per_proc[i] = rcv_buff_size;
			displc[i] = i * rcv_buff_size;
		}
		elements_per_proc[comm_size-1] += numNums % comm_size;
	}

	int rcv_buff[rcv_buff_size];
	MPI_Scatterv(numsTruthSort, elements_per_proc, displc, MPI_INT, rcv_buff,rcv_buff_size,MPI_INT,0,MPI_COMM_WORLD);

	int minMax[2];

	getMinMax(rcv_buff,rcv_buff_size,minMax);

	int global_min, global_max;

	MPI_Reduce(&minMax[0], &global_min, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
	MPI_Reduce(&minMax[1], &global_max, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);

	if (proc_num == 0) {

		sprintf(strBuff0 + strlen(strBuff0),"\nRandomly Generated Numbers:\n");
		printArray(numsTruthSort, numNums, strBuff0);

//		printf("\nNums for sample sort:\n");
//		printArray(inputSS, numNums);

		bubblesort(numsTruthSort,0,numNums);

		time1 = microtime();
		samplesort(inputSS, 0, numNums, numOfProcesses, samplesPerBucket);

		time2 = microtime();

		sprintf(strBuff0 + strlen(strBuff0), "\nSorted Numbers:\n");
		printArray(inputSS, numNums, strBuff0);

		int sum_abs_diff = compareArray(numsTruthSort,inputSS, numNums);
		sprintf(strBuff0 + strlen(strBuff0), "sum of abs value of diff = %d \n\n",sum_abs_diff);
		t = time2 - time1;

		sprintf(strBuff0 + strlen(strBuff0), "\nTime = %g us\n", t);
		sprintf(strBuff0 + strlen(strBuff0), "Timer Resolution = %g us\n", getMicrotimeResolution());

	}

	for (int i = 0; i < comm_size; ++i) {
		if (proc_num == i)
			sprintf(strBuff + strlen(strBuff), "\n proc # %d : elements_per_proc %d\n", proc_num, rcv_buff_size);
	}

	for (int i = 0; i < comm_size; ++i) {
		if (proc_num == i) {
			sprintf(strBuff + strlen(strBuff), "\nRcv buffer:\n");
			printArray(rcv_buff,rcv_buff_size,strBuff);
		}
	}

	sprintf(strBuff + strlen(strBuff), "\nlocal min: %d, local max: %d\n",minMax[0], minMax[1]);

	if (proc_num == 0) {
		sprintf(strBuff0 + strlen(strBuff0), "\nglobal min: %d, global max: %d\n\n",global_min, global_max);
	}

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

	if (proc_num == 0) {
//		printArray(str_buffer_sizes, comm_size, strBuff0);
		printf("\n\n------------------------------------------------------------------\n\n");
		printf("%s\n",strBuff0);
		printf("\n\n------------------------------------------------------------------\n\n");
		printf("%s\n",outputString);
	}

	for (int i = 0; i < comm_size; ++i) {
		if (proc_num == i) {
//			printf("%s",strBuff);
		}
	}
	
	free(outputString);
	MPI_Finalize();
	return 0;
}
