#include "samplesortCPPThreads.hpp"

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

void buildBuckets(std::vector<int>* valsToSort, std::vector<std::vector<int>>* buckets, std::vector<std::vector<int>>* splitters,
    long long offset, long long indexRange, int numOfBuckets, int samplesPerBucket, std::mutex* mtx)
{
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
                    std::lock_guard<std::mutex> lk(*mtx);
                    (*splitters)[i][0] = valToSort;
                    (*buckets)[i].push_back(valToSort);
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
                    std::lock_guard<std::mutex> lk(*mtx);
                    (*splitters)[i][samplesPerBucket - 1] = valToSort;
                    (*buckets)[i].push_back(valToSort);
                    break;
                }
                else
                {
                    wasTooBig = true;
                }
            }
            else
            {
                std::lock_guard<std::mutex> lk(*mtx);
                (*buckets)[i].push_back(valToSort);
                break;
            }
        }
    }
}

void sortAndCombine(std::vector<int>* bucket, std::vector<int>* putHere, long long offset, long long numToDo)
{
    int index = 0;
    std::sort(bucket->begin(), bucket->end());
    for (int i = offset; i < offset + numToDo; i++)
    {
        (*putHere)[i] = (*bucket)[index++];
    }
}

void samplesort(std::vector<int>* toSort, int samplesPerBucket, int numOfBuckets)
{
    std::vector<std::thread*> threads;
    std::mutex mtx;
    int n = toSort->size();
    int numProcessors = numOfBuckets;

    // Get and sort samples
    std::vector<int> samples;
    for (int i = 0; i < samplesPerBucket * numOfBuckets - 1; i++)
    {
        int val = (*toSort)[i * (toSort->size() / (samplesPerBucket * numOfBuckets))];
        //int val = (*toSort)[i * (toSort->size() / (samplesPerBucket * numOfBuckets))];
        samples.push_back(val);
    }
    samples.push_back((*toSort)[toSort->size() - 1]);

    std::sort(samples.begin(), samples.end());



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


    std::vector<std::vector<int>> buckets(numOfBuckets);
    for (long long i = 0; i < numProcessors; i++)
    {
        threads.push_back(new std::thread(buildBuckets, toSort, &buckets, &splitters, static_cast<long long> (i * n / numProcessors),
            static_cast<long long> (n / numProcessors), numOfBuckets, samplesPerBucket, &mtx));
    }

    for (auto t : threads)
    {
        t->join();
    }

    threads.clear();


    /*for (auto valToSort : *toSort)
    {
        bool wasTooBig = false;
        bool wasTooSmall = false;

        for (int i = 0; i < numOfBuckets; i++)
        {
            if (valToSort < splitters[i][0])
            {
                if (wasTooBig || i == 0)
                {
                    splitters[i][0] = valToSort;
                    buckets[i].push_back(valToSort);
                    break;
                }
                else
                {
                    wasTooSmall = true;
                }
            }
            else if (valToSort > splitters[i][samplesPerBucket - 1])
            {
                if (wasTooSmall || i == numOfBuckets - 1)
                {
                    splitters[i][samplesPerBucket - 1] = valToSort;
                    buckets[i].push_back(valToSort);
                    break;
                }
                else
                {
                    wasTooBig = true;
                }
            }
            else
            {
                buckets[i].push_back(valToSort);
                break;
            }
        }
    }*/

    int index = 0;
    int lastBucketSize = 0;
    for (int i = 0; i < numOfBuckets; i++)
    {
        std::vector<int>* bucket = &(buckets[i]);
        int numToDo = bucket->size();
        threads.push_back(new std::thread(sortAndCombine, bucket, toSort, lastBucketSize, numToDo));
        lastBucketSize += buckets[i].size();
    }

    for (auto t : threads)
    {
        t->join();
    }

    threads.clear();
}

void generateVector(std::vector<int>* vec, long long offset, long long indexRange, int randSeed, int upperBound)
{
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

    // Set seed if available
    if (argc >= 5)
    {
        try
        {
            srand(std::stoi(argv[4]));
        }
        catch (...)
        {
            // Value entered is not an int
        }
    }

    std::vector<int> toSort(numNums);

    // Speed up the building of the sortable vector
    std::vector<std::thread*> assignThreads;
    for (long long i = 0; i < numOfProcessors; i++)
    {
        // Perform additional randomization for the vector assignments or else we would just end up with a repeating pattern
        int deepRandSeed = getRandInt(0, 5 * numOfProcessors);
        assignThreads.push_back(new std::thread(generateVector, &toSort, i * numNums / numOfProcessors, numNums / numOfProcessors, deepRandSeed, randomNumBound));
    }
    for (auto t : assignThreads)
    {
        t->join();
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
