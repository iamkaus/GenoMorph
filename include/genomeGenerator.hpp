#pragma once

#include "regionGenerator.hpp"

#include <vector>
#include <string>
#include <random>
#include <ctime>

/**
 * @struct BaseInfo
 * @brief conatins information regard position of base in the sequence as well as base type amongst [A,T,G,C]
 */

struct BaseInfo {

    char        base;
    size_t      position;
};

class GenomeGenerator {
private:
    std::mt19937 rng; /**< Random number generator seeded with system clock. */

    char generate_base(RegionInfo region);

public:

    GenomeGenerator();

    std::vector<BaseInfo> generate_sequence(size_t currentGenomeLength, size_t length);

    std::vector<BaseInfo> complementary_strand(const std::vector<BaseInfo> &original); 

};
