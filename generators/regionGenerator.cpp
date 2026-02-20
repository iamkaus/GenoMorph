#include "regionGenerator.hpp"
#include <random>
#include <fstream>
#include <iostream>
#include <iomanip>

RegionGenerator::RegionGenerator() {
    rng.seed(static_cast<unsigned>(time(0))); // Seed with system clock
}

RegionInfo RegionGenerator::createRegion(size_t currentGenomeLength, size_t genomeLength) {
}

std::array<double, 4> RegionGenerator::regionBasedBaseProbabilities(const RegionInfo &region) {
}