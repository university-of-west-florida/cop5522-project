#include "samplesortPThreads.hpp"

using namespace std;

void printVector(vector<int> vec)
{
    for (long unsigned int i = 0; i < vec.size(); i++)
    {
        cout << to_string(i) << ": " << to_string(vec[i]) << "\n";
    }
    cout << endl;
}

int getRandInt(int min, int max)
{
    int range = max - min;
    if (range == 0) return 0;
    return rand() % range + min;
}

void* buildBuckets(void* vArgs)
{
    buildBucketsArgs* bArgs = (buildBucketsArgs*)vArgs;
    vector<int>* valsToSort             = bArgs->valsToSort;
    vector<vector<int>>* buckets        = bArgs->buckets;
    vector<vector<int>>* splitters      = bArgs->splitters;
    long long offset                    = bArgs->offset;
    long long indexRange                = bArgs->indexRange;
    int numOfBuckets                    = bArgs->numOfBuckets;
    int samplesPerBucket                = bArgs->samplesPerBucket;
    pthread_mutex_t* mtx                = bArgs->mtx;

    for (long long i = offset; i < offset + indexRange; i++)
    {
        int valToSort = (*valsToSort)[i];
        bool wasTooBig = false;
        bool wasTooSmall = false;

        for (int j = 0; j < numOfBuckets; j++)
        {
            if (valToSort < (*splitters)[j][0])
            {
                if (wasTooBig || j == 0)
                {
                    pthread_mutex_lock(mtx);
                    (*splitters)[j][0] = valToSort;
                    (*buckets)[j].push_back(valToSort);
                    pthread_mutex_unlock(mtx);
                    break;
                }
                else
                {
                    wasTooSmall = true;
                }
            }
            else if (valToSort > (*splitters)[j][samplesPerBucket - 1])
            {
                if (wasTooSmall || j == numOfBuckets - 1)
                {
                    pthread_mutex_lock(mtx);
                    (*splitters)[j][samplesPerBucket - 1] = valToSort;
                    (*buckets)[j].push_back(valToSort);
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
                (*buckets)[j].push_back(valToSort);
                pthread_mutex_unlock(mtx);
                break;
            }
        }
    }

    return nullptr;
}

void* sortAndCombine(void* vArgs)//vector<int>* bucket, vector<int>* putHere, long long offset, long long numToDo)
{
    sortAndCombineArgs* sArgs = (sortAndCombineArgs*)vArgs;

    vector<int>* bucket  = sArgs->bucket;
    vector<int>* putHere = sArgs->putHere;
    long long offset          = sArgs->offset;
    long long numToDo         = sArgs->numToDo;

    int index = 0;
    sort(bucket->begin(), bucket->end());
    for (int i = offset; i < offset + numToDo; i++)
    {
        (*putHere)[i] = (*bucket)[index++];
    }

    return nullptr;
}

void samplesort(vector<int>* toSort, int samplesPerBucket, int numOfBuckets)
{
    vector<pthread_t*> threads;
    pthread_mutex_t mtx = PTHREAD_MUTEX_INITIALIZER;
    pthread_mutex_init(&mtx, NULL);
    int n = toSort->size();

    // Get and sort samples
    vector<int> samples;
    for (int i = 0; i < samplesPerBucket * numOfBuckets - 1; i++)
    {
        int val = (*toSort)[i * (toSort->size() / (samplesPerBucket * numOfBuckets))];
        samples.push_back(val);
    }
    samples.push_back((*toSort)[toSort->size() - 1]);

    sort(samples.begin(), samples.end());


    // Get splitters and build out starting buckets
    vector<vector<int>> splitters;
    for (int i = 0; i < numOfBuckets; i++)
    {
        vector<int> split;
        for (int j = 0; j < samplesPerBucket; j++)
        {
            split.push_back(samples[i * samplesPerBucket + j]);
        }
        splitters.push_back(split);
    }

    // Build out our buckets with our splitters
    vector<vector<int>> buckets(numOfBuckets);
    for (long long i = 0; i < 1; i++)
    {
        buildBucketsArgs* bArgs = new buildBucketsArgs();
        bArgs->valsToSort       = toSort;
        bArgs->buckets          = &buckets;
        bArgs->splitters        = &splitters;
        bArgs->offset           = static_cast<long long> (i * n / 1);
        bArgs->indexRange       = static_cast<long long> (n / 1);
        bArgs->numOfBuckets     = numOfBuckets;
        bArgs->samplesPerBucket = samplesPerBucket;
        bArgs->mtx              = &mtx;

        buildBuckets(static_cast<void*>(bArgs));
        /*pthread_t* newThread = new pthread_t();
        threads.push_back(newThread);

        if (pthread_create(newThread, NULL, buildBuckets, static_cast<void*>(bArgs)))
        {
            cout << "Pthread_create error: " << strerror(errno) << endl;
            exit(1);
        }*/
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
    int lastBucketSize = 0;
    for (int i = 0; i < numOfBuckets; i++)
    {
        vector<int>* bucket = &(buckets[i]);
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
            cout << "Pthread_create error: " << strerror(errno) << endl;
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
    vector<int>* vec = gArgs->vec;
    long long offset     = gArgs->offset;
    long long indexRange = gArgs->indexRange;
    int randSeed         = gArgs->randSeed;
    int upperBound       = gArgs->upperBound;

    srand(randSeed);
    for (long long i = offset; i < offset + indexRange; i++)
    {
        (*vec)[i] = getRandInt(0, upperBound);
    }

    return nullptr;
}


int main(int argc, char** argv)
{
    if (argc < 4) {
        fprintf(stderr, "USAGE: %s <numNums> <randomNumBound> <samplesPerBucket> (opt)<random seed (-1 for none)> (opt)<number of processors>\n", argv[0]);
        exit(1);
    }

    double t, time1, time2;

    // Parsing arguments
    int numNums = atoi(argv[1]);
    int randomNumBound = atoi(argv[2]);
    int samplesPerBucket = atoi(argv[3]);
    int numOfProcessors;


    if (argc >= 6)
    {
        numOfProcessors = stoi(argv[5]);
    }
    else
    {
        // Getting the current hardware's number of logical processors
        numOfProcessors = sysconf(_SC_NPROCESSORS_ONLN);
    }
    cout << "Numofprocessors = " << to_string(numOfProcessors) << endl;
    

    vector<int> toSort(numNums);

    // Set seed if available
    if (argc >= 5)
    {
        if (string(argv[4]) != "-1")
        {
            try
            {
                srand(stoi(argv[4]));
                // Serially add in random numbers to the vector
                for (int i = 0; i < numNums; i++)
                    toSort[i] = getRandInt(0, randomNumBound);
            }
            catch (...)
            {
                // Value entered is not an int -- maybe alert user?
                cout << "ERROR: Entered invalid seed" << endl;
                exit(1);
            }
        }
    }
    // Else, use threads to speed up the building of the vector with extra random elements
    else
    {
        // Speed up the building of the sortable vector
        vector<pthread_t*> assignThreads;
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
                cout << "Pthread_create error: " << strerror(errno) << endl;
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

    for (unsigned long int i = 1; i < toSort.size(); i++)
    {
        if (toSort[i - 1] > toSort[i])
        {
            success = false;
            cout << "Failed at i: " << to_string(i) << endl;
        }
    }

    string result = (success) ? "PASSED" : "FAILED";
    cout << result << endl;

    return 0;
}
