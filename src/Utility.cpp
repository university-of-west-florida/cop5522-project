#include "Utility.hpp"

// Print any type of vector
template<typename T>
void printVector(std::vector<T>& toPrint)
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


template<typename T>
void printAny(T toPrint, std::string name)
{
    try
    {
        std::cout << name << " = " << std::to_string(toPrint) << "\t";
    }
    catch (...)
    {
        std::cout << "ERROR: Unable to print variable: " << name << std::endl;
    }
}

void printDebug()
{
    std::cout << std::endl;
}

template<typename T, typename... Args>
void printDebug(T toPrint, Args... args)
{
    std::string name = currentDebugNames.substr(0, currentDebugNames.find(","));
    currentDebugNames.erase(0, currentDebugNames.find(",") + 1);
    printAny(toPrint, name);
    printDebug(args...);
};
