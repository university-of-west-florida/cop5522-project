#include "samplesortPThreads.hpp"

using namespace std;

void printVector(std::vector<int> vec)
{
    for (int i = 0; i < vec.size(); i++)
    {
        std::cout << std::to_string(i) << ": " << std::to_string(vec[i]) << "\n";
    }
    std::cout << std::endl;
}

int getRandInt(int min, int max)
{
    int range = max - min;
    if (range == 0) return 0;
    return rand() % range + min;
}

int ceiling(int x, int y)
{
    return (x + y - 1) / y;
}

void* buildBuckets(void* vArgs)
{
    buildBucketsArgs* bArgs = (buildBucketsArgs*)vArgs;
    std::vector<int>* valsToSort             = bArgs->valsToSort;
    std::vector<std::vector<int>>* buckets   = bArgs->buckets;
    std::vector<std::vector<int>>* splitters = bArgs->splitters;
    long long offset                         = bArgs->offset;
    long long indexRange                     = bArgs->indexRange;
    int numOfBuckets                         = bArgs->numOfBuckets;
    int samplesPerBucket                     = bArgs->samplesPerBucket;
    pthread_mutex_t* mtx                     = bArgs->mtx;

    for (long long i = offset; i < offset + indexRange; i++)
    {
        int valToSort = (*valsToSort)[i];
        bool wasTooBig = false;
        bool wasTooSmall = false;

        for (int i = 0; i < numOfBuckets; i++)
        {
            if (valToSort < (*splitters)[i][0])
            {
                if (wasTooBig || i == 0)
                {
                    pthread_mutex_lock(mtx);
                    (*splitters)[i][0] = valToSort;
                    (*buckets)[i].push_back(valToSort);
                    pthread_mutex_unlock(mtx);
                    break;
                }
                else
                {
                    wasTooSmall = true;
                }
            }
            else if (valToSort > (*splitters)[i][samplesPerBucket - 1])
            {
                if (wasTooSmall || i == numOfBuckets - 1)
                {
                    pthread_mutex_lock(mtx);
                    (*splitters)[i][samplesPerBucket - 1] = valToSort;
                    (*buckets)[i].push_back(valToSort);
                    pthread_mutex_unlock(mtx);
                    break;
                }
                else
                {
                    wasTooBig = true;
                }
            }
            else
            {
                pthread_mutex_lock(mtx);
                (*buckets)[i].push_back(valToSort);
                pthread_mutex_unlock(mtx);
                break;
            }
        }
    }
}

void* sortAndCombine(void* vArgs)//std::vector<int>* bucket, std::vector<int>* putHere, long long offset, long long numToDo)
{
    sortAndCombineArgs* sArgs = (sortAndCombineArgs*)vArgs;

    std::vector<int>* bucket  = sArgs->bucket;
    std::vector<int>* putHere = sArgs->putHere;
    long long offset          = sArgs->offset;
    long long numToDo         = sArgs->numToDo;

    int index = 0;
    std::sort(bucket->begin(), bucket->end());
    for (int i = offset; i < offset + numToDo; i++)
    {
        (*putHere)[i] = (*bucket)[index++];
    }
}

void samplesort(std::vector<int>* toSort, int samplesPerBucket, int numOfBuckets)
{
    std::vector<pthread_t*> threads;
    pthread_mutex_t mtx = PTHREAD_MUTEX_INITIALIZER;
    pthread_mutex_init(&mtx, NULL);
    int n = toSort->size();
    int numProcessors = numOfBuckets;

    // Get and sort samples
    std::vector<int> samples;
    for (int i = 0; i < samplesPerBucket * numOfBuckets - 1; i++)
    {
        int val = (*toSort)[i * (toSort->size() / (samplesPerBucket * numOfBuckets))];
        samples.push_back(val);
    }
    samples.push_back((*toSort)[toSort->size() - 1]);

    std::sort(samples.begin(), samples.end());


    // Get splitters and build out starting buckets
    std::vector<std::vector<int>> splitters;
    for (int i = 0; i < numOfBuckets; i++)
    {
        std::vector<int> split;
        for (int j = 0; j < samplesPerBucket; j++)
        {
            split.push_back(samples[i * samplesPerBucket + j]);
        }
        splitters.push_back(split);
    }

    // Build out our buckets with our splitters
    std::vector<std::vector<int>> buckets(numOfBuckets);
    for (long long i = 0; i < numProcessors; i++)
    {
        buildBucketsArgs* bArgs = new buildBucketsArgs();
        bArgs->valsToSort       = toSort;
        bArgs->buckets          = &buckets;
        bArgs->splitters        = &splitters;
        bArgs->offset           = static_cast<long long> (i * n / numProcessors);
        bArgs->indexRange       = static_cast<long long> (n / numProcessors);
        bArgs->numOfBuckets     = numOfBuckets;
        bArgs->samplesPerBucket = samplesPerBucket;
        bArgs->mtx              = &mtx;

        pthread_t* newThread = new pthread_t();
        threads.push_back(newThread);

        if (pthread_create(newThread, NULL, buildBuckets, static_cast<void*>(bArgs)))
        {
            std::cout << "Pthread_create error: " << strerror(errno) << std::endl;
            exit(1);
        }
    }

    // Wait for all threads to finish building out the buckets
    // and delete and clear them from the vector
    for (auto t : threads)
    {
        pthread_join(*t, NULL);
        delete[] t;
    }
    threads.clear();


    // Sort and combine our buckets back to original vector
    int index = 0;
    int lastBucketSize = 0;
    for (int i = 0; i < numOfBuckets; i++)
    {
        std::vector<int>* bucket = &(buckets[i]);
        int numToDo = bucket->size();

        sortAndCombineArgs* sArgs = new sortAndCombineArgs();
        sArgs->bucket  = bucket;
        sArgs->putHere = toSort;
        sArgs->offset  = lastBucketSize;
        sArgs->numToDo = numToDo;

        pthread_t* newThread = new pthread_t();
        threads.push_back(newThread);

        if (pthread_create(newThread, NULL, sortAndCombine, static_cast<void*>(sArgs)))
        {
            std::cout << "Pthread_create error: " << strerror(errno) << std::endl;
            exit(1);
        }
        lastBucketSize += buckets[i].size();
    }

    // Wait for all threads to finish their sorting and combining back to original vector
    // and delete and clear them from the vector
    for (const auto& t : threads)
    {
        pthread_join(*t, NULL);
        delete[] t;
    }
    threads.clear();
}

void* generateVector(void* vArgs)
{
    generateVectorArgs* gArgs = (generateVectorArgs*)vArgs;
    std:vector<int>* vec = gArgs->vec;
    long long offset     = gArgs->offset;
    long long indexRange = gArgs->indexRange;
    int randSeed         = gArgs->randSeed;
    int upperBound       = gArgs->upperBound;

    srand(randSeed);
    for (long long i = offset; i < offset + indexRange; i++)
    {
        (*vec)[i] = getRandInt(0, upperBound);
    }
}


int main(int argc, char** argv)
{
    if (argc < 4) {
        fprintf(stderr, "USAGE: %s <numNums> <randomNumBound> <samplesPerBucket> (opt)<random seed>\n", argv[0]);
        exit(1);
    }

    double t, time1, time2;

    // Parsing arguments
    int numNums = atoi(argv[1]);
    int randomNumBound = atoi(argv[2]);
    int samplesPerBucket = atoi(argv[3]);

    // Getting the current hardware's number of logical processors
    int numOfProcessors = sysconf(_SC_NPROCESSORS_ONLN);
    

    std::vector<int> toSort(numNums);

    // Set seed if available
    if (argc >= 5)
    {
        try
        {
            srand(std::stoi(argv[4]));
            // Serially add in random numbers to the vector
            for (int i = 0; i < numNums; i++)
                toSort[i] = getRandInt(0, randomNumBound);
        }
        catch (...)
        {
            // Value entered is not an int -- maybe alert user?
            std::cout << "ERROR: Entered invalid seed" << std::endl;
            exit(1);
        }
    }
    // Else, use threads to speed up the building of the vector with extra random elements
    else
    {
        // Speed up the building of the sortable vector
        std::vector<pthread_t*> assignThreads;
        for (long long i = 0; i < numOfProcessors; i++)
        {
            // Perform additional randomization for the vector assignments or else we would just end up with a repeating pattern
            int deepRandSeed = getRandInt(0, 5 * numOfProcessors);

            generateVectorArgs* gArgs = new generateVectorArgs();
            gArgs->vec = &toSort;
            gArgs->offset = i * numNums / numOfProcessors;
            gArgs->indexRange = numNums / numOfProcessors;
            gArgs->randSeed = deepRandSeed;
            gArgs->upperBound = randomNumBound;

            pthread_t* newThread = new pthread_t();
            assignThreads.push_back(newThread);

            if (pthread_create(newThread, NULL, generateVector, static_cast<void*>(gArgs)))
            {
                std::cout << "Pthread_create error: " << strerror(errno) << std::endl;
                exit(1);
            }
        }
        for (auto t : assignThreads)
        {
            pthread_join(*t, NULL);
            delete[] t;
        }
        assignThreads.clear();
    }

    // measure time
    time1 = microtime();
    samplesort(&toSort, samplesPerBucket, numOfProcessors);
    time2 = microtime();
    t = time2 - time1;


    // Print results
    printf("\nTime = %g us\n", t);
    printf("Timer Resolution = %g us\n", getMicrotimeResolution());


    // Determine if sort was successful
    bool success = true;

    for (int i = 1; i < toSort.size(); i++)
    {
        if (toSort[i - 1] > toSort[i])
        {
            success = false;
            std::cout << "Failed at i: " << std::to_string(i) << std::endl;
        }
    }

    std::string result = (success) ? "PASSED" : "FAILED";
    std::cout << result << std::endl;

    return 0;
}
