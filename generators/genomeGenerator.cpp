#include "genomeGenerator.hpp"
#include <iostream>
#include <sstream>
#include <random>
#include <ctime>
#include <fstream>


GenomeGenerator::GenomeGenerator()
{
    rng.seed(static_cast<unsigned>(time(0))); // Seed with system clock
}

char GenomeGenerator::generate_base(RegionInfo region) {
}

std::vector<BaseInfo> GenomeGenerator::generate_sequence(size_t total_generated, size_t length) {
}