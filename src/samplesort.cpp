#include <microtime.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <vector>
#include <time.h>

void printArray(int inputNums[], int numNums) {
    printf("Array Values:\n\t");
    int i;
    for(i = 0; i < numNums; ++i)
        printf ("%d ",inputNums[i]);
    printf("\n\n");
}

// Initializing array with numNums values which are below randomNumBound
void generateRandNums(int inputNums[], int numNums, int randomNumBound) {
    // Init random number generator the assign them with the specified bound
    int i;
    srand(time(0));
    for(i = 0; i < numNums; ++i)
        inputNums[i] = rand() % randomNumBound;
}

// https://stackoverflow.com/questions/17655748/bubble-sort-algorithm-in-c
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

void samplesort(int inputNums[], int numNums, int numOfBuckets, int samplesPerBucket) {
    int i, j, k;

    // If there are not enough inputNums to take the quantity of samples per bucket, then just set samplesPerBucket to 1
    if(numOfBuckets * samplesPerBucket > numNums) {
        if(numNums / numOfBuckets >= 1) {
            samplesPerBucket = 1;
        }
        // If there are fewer inputNums than specified numOfBuckets, then just bubblesort the rest bc it is likely a small set left
        else {
            bubblesort(inputNums, 0, numNums);
            return;
        }
    }

    // Creating an array of samples so that we can find pivot points for the algorithm
    int sortedSamples[(numOfBuckets - 1) * samplesPerBucket];
    for(i = 0; i < (numOfBuckets - 1) * samplesPerBucket; ++i)
        sortedSamples[i] = inputNums[i];

    // Sorting the collection of samples with bubble sort
    bubblesort(sortedSamples, 0, (numOfBuckets - 1) * samplesPerBucket);
    
    // Condensing the sortedSamples into just splitters between the buckets
    int splitters[numOfBuckets - 1];
    for(i = 0; i < numOfBuckets - 1; ++i)
        splitters[i] = sortedSamples[i * samplesPerBucket];

    // Going through inputNums and sorting each value into a bucket (each bucket is stored as a vector)
    // NOTE: The initial double nesting of for loops is not necessary for serial purposes. I could have made that a single
    //       for loop, but it will be easier to parallelize it when set up this way
    std::vector<int> buckets[numOfBuckets];
    for(i = 0; i < numOfBuckets; ++i) {
        for(j = i * (numNums / numOfBuckets); j < (i + 1) * (numNums / numOfBuckets); ++j) {
            // Looping through to check which bucket the current value belongs to, 
            // then pushing it into the bucket if it falls in that bound
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

    // Handling the remainder values at the end of the list which were not accounted for above
    for(i = numNums - (numNums % numOfBuckets); i < numNums; ++i) {
        for(k = 0; k < numOfBuckets - 1; ++k) {
            if(inputNums[j] <= splitters[k]) {
                buckets[k].push_back(inputNums[j]);
                break;
            }
        }
        buckets[numOfBuckets - 1].push_back(inputNums[j]);
    }

    // Goes one-by-one through each bucket and copies their values into inputNums
    int runningIndex = 0;
    for(i = 0; i < numOfBuckets; ++i) {
        int bucketSize = buckets[i].size();
        for(j = 0; j < bucketSize; ++j ) {
            inputNums[runningIndex + j] = buckets[i].at(j);
        }
        // Once the values are copied into inputNums, each bucket is sorted in place (using bubble sort for now)
        bubblesort(inputNums, runningIndex, runningIndex + bucketSize);
        runningIndex += bucketSize;
    }
}

int main(int argc, char** argv) {
    double t, time1, time2;

    if (argc != 4) {
        fprintf(stderr, "USAGE: %s <numNums> <randomNumBound> <samplesPerBucket>\n", argv[0]);
        exit(1);
    }

    int numNums = atoi(argv[1]);
    int randomNumBound = atoi(argv[2]);
    int samplesPerBucket = atoi(argv[3]);

    int inputNums[numNums];
    generateRandNums(inputNums, numNums, randomNumBound);
    printArray(inputNums, numNums);

    int numOfProcessors = sysconf(_SC_NPROCESSORS_ONLN);

    // measure time
    time1 = microtime();
    samplesort(inputNums, numNums, numOfProcessors, samplesPerBucket);
    //bubblesort(inputNums, 0, numNums);
    time2 = microtime();

    printArray(inputNums, numNums);

    t = time2 - time1;

    // Print results
    printf("\nTime = %g us\n", t);
    printf("Timer Resolution = %g us\n", getMicrotimeResolution());

    return 0;
}