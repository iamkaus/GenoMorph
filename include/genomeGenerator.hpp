#pragma once

#include "regionGenerator.hpp"

#include <vector>
#include <string>
#include <random>
#include <ctime>

struct BaseInfo {
 
    char        base;
    size_t      position;
}; // maybe inherit from RegionInfo struct from "regionGenerator.hpp".

class GenomeGenerator {
private:
    std::mt19937 rng; /**< Random number generator seeded with system clock. */

    char generate_base(RegionInfo region);

public:

    GenomeGenerator();

    std::vector<BaseInfo> generate_sequence(size_t currentGenomeLength, size_t length);

    std::vector<BaseInfo> complementary_strand(const std::vector<BaseInfo> &original); 

};
