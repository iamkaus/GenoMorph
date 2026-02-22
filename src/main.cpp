#include "genomeGenerator.hpp"
#include "regionGenerator.hpp"

int main() {

    GenomeGenerator generator;

    auto sequence = generator.generate_sequence(0, 10000);

    return 0;
}