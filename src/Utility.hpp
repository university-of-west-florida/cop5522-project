#pragma once

#include <string>

#define DEBUG(...) currentDebugNames = #__VA_ARGS__; printDebug(__VA_ARGS__)
std::string currentDebugNames;
#define ENDL std::cout << std::endl
