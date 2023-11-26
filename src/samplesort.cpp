#include "samplesort.h"


void printArray(int inputNums[], int numNums) {
    printf("Array Values:\n\t");
    int i;
    for(i = 0; i < numNums; ++i)
        printf ("%d ",inputNums[i]);
    printf("\n\n");
}

// Print any type of vector
template<typename T>
void printVector(std::vector<T> &toPrint)
{
    int index = 0;
    for (const auto& it : toPrint)
    {
        try
        {
            std::cout << index++ << ": " << std::to_string(it) << std::endl;
        }
        catch (...)
        {
            std::cout << "ERROR: Unable to print variable. Variable type: " << typeid(toPrint).name() << std::endl;
            return;
        }
    }   
}


// Initializing vector with numNums values which are below randomNumBound
void generateRandNums(std::vector<int> &inputNums, int numNums, int randomNumBound) {
    // Init random number generator the assign them with the specified bound
    int i;
    srand(time(0));
    for (i = 0; i < numNums; ++i)
        inputNums.push_back(rand() % randomNumBound);
}

// Leave original array version for backwards compatibility
// Initializing array with numNums values which are below randomNumBound
void generateRandNums(int inputNums[], int numNums, int randomNumBound) {
    // Init random number generator the assign them with the specified bound
    int i;
    srand(time(0));
    for (i = 0; i < numNums; ++i)
        inputNums[i] = (rand() % randomNumBound);
}

void bubblesort(std::vector<int> &inputNums, int start, int end) {
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

void* testFunc(void* vArgs)
{
    // Convert void* back to struct*
    testArgs* tArgs = static_cast<testArgs*>(vArgs);
    int j = tArgs->j;
    int numOfBuckets = tArgs->numOfBuckets;

    // Copy out vectors for easier access
    std::vector<int>* inputNums = tArgs->inputNums;
    std::vector<int>* splitters = tArgs->splitters;
    std::vector<std::vector<int>>* buckets = tArgs->buckets;

    for (int k = 0; k < numOfBuckets - 1; ++k) {
        if ((*inputNums)[j] <= (*splitters)[k]) {
            (*buckets)[k].push_back((*inputNums)[j]);
            break;
        }
        if (k == numOfBuckets - 2)
        {
            (*buckets)[numOfBuckets - 1].push_back((*inputNums)[j]);
        }
    }
}

void* samplesort(void* vArgs) {
    sampleSortArgs* args = (sampleSortArgs*)vArgs;
    
    std::vector<int>* inputNums = args->inputNums;
    int start = args->start;
    int end = args->end;
    int numOfBuckets = args->numOfBuckets;
    int samplesPerBucket = args->samplesPerBucket;
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
            bubblesort((*inputNums), start, end);
            return nullptr;
        }
    }

    // STEP 2
    // Creating an array of samples so that we can find pivot points for the algorithm
    std::vector<int> sortedSamples((numOfBuckets - 1) * samplesPerBucket);
    for (i = start; i < start + (numOfBuckets - 1) * samplesPerBucket; ++i)
        sortedSamples[i - start] = (*inputNums)[i];

    // Sorting the collection of samples with bubble sort
    bubblesort(sortedSamples, 0, (numOfBuckets - 1) * samplesPerBucket);

    // Condensing the sortedSamples into just splitters between the buckets
    std::vector<int> splitters(numOfBuckets - 1);
    for (i = 0; i < numOfBuckets - 1; ++i)
        splitters[i] = sortedSamples[i * samplesPerBucket];

    // STEP 3
    // Going through inputNums and sorting each value into a bucket (each bucket is stored as a vector)
    // NOTE: The initial double nesting of for loops is not necessary for serial purposes. I could have made the first
    //        two for loops into a single for loop, but I think it will be easier to parallelize it when set up this way
    std::vector<std::vector<int>> buckets(numOfBuckets);
    testArgs tArgs = { inputNums, &splitters, &buckets };
    tArgs.numOfBuckets = numOfBuckets;
    // This breaks the initial inputNums into equal segments, then loops through each value and puts it into the bucket it belongs to
    for (i = 0; i < numOfBuckets; ++i)
    {
        for (j = start + (i * (numNums / numOfBuckets)); j < start + ((i + 1) * (numNums / numOfBuckets)); ++j)
        {
            tArgs.j = j;

            pthread_mutex_lock(&(args->mtx));
            pthread_t* newThread = new pthread_t();

            // Thread number is for debug only
            // TODO:: Remove after testing
            tArgs.threadNum = numThreads++;

            // Save off threads to join them after being spawned
            args->threads->push_back(newThread);

            // Start up thread
            pthread_create(newThread, NULL, testFunc, static_cast<void*>(& tArgs));
            pthread_mutex_unlock(&(args->mtx));
        }
    }
    // Wait for all of our threads to complete
    for (const auto& thread : *(args->threads))
    {
        pthread_join(*thread, NULL);
    }

    // Handling the remainder values at the end of the list which were not accounted for above
    // The inner for loop is the same as above bc it just sorts the value in question into the right bucket
    // If you have 20 buckets, but 23 values in the array, the above loop will handle the first 20 values
    // and this loop will handle the remaining three bc they are kinda outliers
    for (i = start + numNums - (numNums % numOfBuckets); i < end; ++i) {
        for (k = 0; k < numOfBuckets - 1; ++k) {
            if ((*inputNums)[i] <= splitters[k]) {
                buckets[k].push_back((*inputNums)[i]);
                break;
            }
            if (k == numOfBuckets - 2)
                buckets[numOfBuckets - 1].push_back((*inputNums)[i]);
        }
    }

    // STEP 4
    // Goes one-by-one through each bucket and copies their values into inputNums
    int runningIndex = start;
    for (i = 0; i < numOfBuckets; ++i) {
        int bucketSize = buckets[i].size();
        for (j = 0; j < bucketSize; ++j)
        {
            (*inputNums)[runningIndex + j] = buckets[i].at(j);
        }
        // Once the values of a bucket are copied into inputNums, we
        // recursively call samplesort to sort those values in place.
        // Might remove this after testing if I cannot get the whole algorithm threading working
        sampleSortArgs ssa = { inputNums, args->threads, runningIndex, runningIndex + bucketSize, numOfBuckets, samplesPerBucket, args->mtx };
        samplesort((void*) &ssa);
        
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
    std::vector<pthread_t*> threads;
    pthread_mutex_t mtx;

    // Generating an array of numNums random numbers which are bounded by randomNumBound
    std::vector<int> inputNums;
    generateRandNums(inputNums, numNums, randomNumBound);
    printVector(inputNums);

    // Getting the current hardware's number of logical processors
    int numOfProcessors = sysconf(_SC_NPROCESSORS_ONLN);

    // measure time
    time1 = microtime();
    sampleSortArgs ssa = { &inputNums, &threads, 0, numNums, numOfProcessors, samplesPerBucket, mtx };
    samplesort((void*) &ssa);
    //samplesort(inputNums, 0, numNums, numOfProcessors, samplesPerBucket);
    // Try out the below bubblesort instead if you want to compare times
    //bubblesort(inputNums, 0, numNums);
    time2 = microtime();

    //printArray(inputNums, numNums);
    printVector(inputNums);

    t = time2 - time1;

    // Print results
    printf("\nTime = %g us\n", t);
    printf("Timer Resolution = %g us\n", getMicrotimeResolution());

    return 0;
}
