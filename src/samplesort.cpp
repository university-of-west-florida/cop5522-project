#include "samplesort.h"

void printArray(int inputNums[], int numNums) {
    printf("Array Values:\n\t");
    int i;
    for (i = 0; i < numNums; ++i)
        printf("%d ", inputNums[i]);
    printf("\n\n");
}

// Initializing array with numNums values which are below randomNumBound
void generateRandNums(int inputNums[], int numNums, int randomNumBound) {
    // Init random number generator the assign them with the specified bound
    int i;
    srand(time(0));
    for (i = 0; i < numNums; ++i)
        inputNums[i] = rand() % randomNumBound;
}

// https://stackoverflow.com/questions/17655748/bubble-sort-algorithm-in-c
// This can be used to compare to if you want. Samplesort is MUCH faster on big sorts
void bubblesort(int inputNums[], int start, int end) {
    int i, j, temp;
    for (i = start; i < end - 1; ++i) {
        for (j = start; j < start + end - i - 1; ++j) {
            if (inputNums[j] > inputNums[j + 1]) {
                temp = inputNums[j + 1];
                inputNums[j + 1] = inputNums[j];
                inputNums[j] = temp;
            }
        }
    }
}

// Baseline Samplesort
//
// STEP 1: Since this is recursive, bubblesort does the sorting at the innermost call when numNums is less than numOfBuckets
//         Drop the number of samplesPerBucket to 1 if that allows the samplesort to continue. Otherwise, just call bubblesort on the rest
//
// STEP 2: Take numOfBuckets * samplesPerBucket samples from inputNums (since inputNums is assumed to be random, I just pull values from the front)
//         Sort these samples using bubblesort (could use samplesort here but it is likely a small array anyway so just bubblesort it)
//         Take numOfBuckets - 1 values evenly from the samples. These are your splitters and should vaguely divide the 
//         entire inputNums into numOfBuckets buckets
//
// STEP 3: Sort all values in inputNums into a bucket (a bucket is defined by two adjacent splitters)
//         For parallel implementation, it can be beneficial to segment inputNums evenly up and stick a processor on each one
//         This is why the implementation segments the values numOfBuckets ways (bc numOfProcessors is what I pass into the numOfBuckets argument)
//
// STEP 4: Copy all buckets into inputNums in order. Sort each bucket in inputNums
//         Once each bucket is sorted in place, there is no work to do to combine the buckets because they are already in order
//         Samplesort is recursively called to sort each bucket in place

void samplesort(int inputNums[], int start, int end, int numOfBuckets, int samplesPerBucket) {
    int i, j, k;
    // numNums is the number of numbers that we are sorting
    int numNums = end - start;

    // STEP 1
    // If there are not enough inputNums to take the quantity of samples per bucket, then just set samplesPerBucket to 1
    if (numOfBuckets * samplesPerBucket > numNums) {
        if (numNums / numOfBuckets >= 1) {
            samplesPerBucket = 1;
        }
        // If there are fewer inputNums than specified numOfBuckets, then just bubblesort the rest bc it is likely a small set left
        else {
            bubblesort(inputNums, start, end);
            return;
        }
    }

    // STEP 2
    // Creating an array of samples so that we can find pivot points for the algorithm
    int sortedSamples[(numOfBuckets - 1) * samplesPerBucket];
    for (i = start; i < start + (numOfBuckets - 1) * samplesPerBucket; ++i)
        sortedSamples[i - start] = inputNums[i];

    // Sorting the collection of samples with bubble sort
    bubblesort(sortedSamples, 0, (numOfBuckets - 1) * samplesPerBucket);

    // Condensing the sortedSamples into just splitters between the buckets
    int splitters[numOfBuckets - 1];
    for (i = 0; i < numOfBuckets - 1; ++i)
        splitters[i] = sortedSamples[i * samplesPerBucket];

    // STEP 3
    // Going through inputNums and sorting each value into a bucket (each bucket is stored as a vector)
    // NOTE: The initial double nesting of for loops is not necessary for serial purposes. I could have made the first
    //        two for loops into a single for loop, but I think it will be easier to parallelize it when set up this way
    std::vector<int> buckets[numOfBuckets];
    // This breaks the initial inputNums into equal segments, then loops through each value and puts it into the bucket it belongs to
    for (i = 0; i < numOfBuckets; ++i) {
        for (j = start + (i * (numNums / numOfBuckets)); j < start + ((i + 1) * (numNums / numOfBuckets)); ++j) {
            // Looping through to check which bucket the current value belongs to, 
            // then pushing it into the bucket if it falls in that bound
            for (k = 0; k < numOfBuckets - 1; ++k) {
                if (inputNums[j] <= splitters[k]) {
                    buckets[k].push_back(inputNums[j]);
                    break;
                }
                if (k == numOfBuckets - 2)
                    buckets[numOfBuckets - 1].push_back(inputNums[j]);
            }
        }
    }

    // Handling the remainder values at the end of the list which were not accounted for above
    // The inner for loop is the same as above bc it just sorts the value in question into the right bucket
    // If you have 20 buckets, but 23 values in the array, the above loop will handle the first 20 values
    // and this loop will handle the remaining three bc they are kinda outliers
    for (i = start + numNums - (numNums % numOfBuckets); i < end; ++i) {
        for (k = 0; k < numOfBuckets - 1; ++k) {
            if (inputNums[i] <= splitters[k]) {
                buckets[k].push_back(inputNums[i]);
                break;
            }
            if (k == numOfBuckets - 2)
                buckets[numOfBuckets - 1].push_back(inputNums[i]);
        }
    }

    // STEP 4
    // Goes one-by-one through each bucket and copies their values into inputNums
    int runningIndex = start;
    for (i = 0; i < numOfBuckets; ++i) {
        int bucketSize = buckets[i].size();
        for (j = 0; j < bucketSize; ++j)
            inputNums[runningIndex + j] = buckets[i].at(j);
        // Once the values of a bucket are copied into inputNums, we
        // recursively call samplesort to sort those values in place.
        samplesort(inputNums, runningIndex, runningIndex + bucketSize, numOfBuckets, samplesPerBucket);
        runningIndex += bucketSize;
    }
}

int main(int argc, char** argv) {
    double t, time1, time2;

    if (argc != 4) {
        fprintf(stderr, "USAGE: %s <numNums> <randomNumBound> <samplesPerBucket>\n", argv[0]);
        exit(1);
    }

    // Parsing arguments
    int numNums = atoi(argv[1]);
    int randomNumBound = atoi(argv[2]);
    int samplesPerBucket = atoi(argv[3]);

    // Generating an array of numNums random numbers which are bounded by randomNumBound
    int inputNums[numNums];
    generateRandNums(inputNums, numNums, randomNumBound);
    //printArray(inputNums, numNums);

    // Getting the current hardware's number of logical processors
    int numOfProcessors = sysconf(_SC_NPROCESSORS_ONLN);

    // measure time
    time1 = microtime();
    samplesort(inputNums, 0, numNums, numOfProcessors, samplesPerBucket);
    // Try out the below bubblesort instead if you want to compare times
    //bubblesort(inputNums, 0, numNums);
    time2 = microtime();

    //printArray(inputNums, numNums);

    t = time2 - time1;

    // Print results
    printf("\nTime = %g us\n", t);
    printf("Timer Resolution = %g us\n", getMicrotimeResolution());

    return 0;
}