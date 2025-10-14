#ifndef REGION_GENERATOR_H
#define REGION_GENERATOR_H

#include <string>
#include <random>
#include <array>

/**
 * @struct RegionState
 * @brief Represents the state and properties of a genomic region.
 * 
 * This struct is used to define regions of DNA with specific characteristics
 * such as type (coding or non-coding), GC content, and length.
*/

struct RegionState {

    /**
     * @brief The starting index of the region in the genome.
     * Zero-based index.
    */
    size_t start_index;

    /**
     * @brief The ending index of the region in the genome.
     * Zero-based index.
    */
    size_t end_index;

    /**
     * @brief Target GC or AT content for the region.
     * Value between 0.35 to 0.45 for GC content and 0.65 to 0.55 for AT content representing the fraction of G and C bases.
    */
    double target_content;

    /**
     * @brief The current GC or AT content generated in the region.
     * Used to track percentage progress towards the target_content.
    */
    double current_content;

    /**
     * @brief The current length of the region to be generated.
     * Number of bases in the region.
    */
    size_t current_length_generated;

    /**
     * @brief The target length of the region to be generated.
     * Number of bases in the region.
    */
    size_t target_length;

    /**
     * @brief The type of the region: "coding" or "non_coding".
     * Determines the characteristics of the region.
    */
    std::string region_type;
};

class RegionGenerator {
private:

    /**
     * @brief Random number generator seeded with system clock.
     * Used for stochastic decisions in region generation.
    */
    std::mt19937 rng;
public:

    /**
     * @brief Constructs a RegionGenerator and seeds the RNG.
     * The RNG is seeded using the current system time to ensure different sequences on each run.
    */
    RegionGenerator();

    /**
     * @brief Creates a new genomic region based on the current genome length and total genome length.
     * @param currentGenomeLength The length of the genome generated so far.
     * @param genomeLength The total desired length of the genome.
     * @return A RegionState struct representing the newly created region.
     * Throws std::invalid_argument if genomeLength is less than 100.
     * The function randomly decides the type of region (coding or non-coding),
     * its target GC or AT content, and its length based on predefined probabilities and distributions.
    */
    RegionState createRegion(size_t currentGenomeLength, size_t genomeLength);

    /**
     * @brief Provides base probabilities based on the region type.
     * @param region The RegionState providing context for base probability determination.
     * @return An array of doubles representing the probabilities for A, T, C, and G respectively.
     * This function adjusts base probabilities to reflect the characteristics of the region.
     * For coding regions, higher GC content is favored, while non-coding regions favor AT content.
     * The returned probabilities can be used in base generation to ensure the sequence adheres to the region's properties.
     * The sum of the returned probabilities should ideally be 1.0.
     * This function is crucial for generating biologically relevant sequences.
     * The design allows for easy modification of base probabilities in the future.
     * Overall, this function encapsulates the logic for region-specific base probability determination.
     */
    std::array<double, 4> regionBasedBaseProbabilities(const RegionState &region);

    /**
     * --------------------------------------------------------------------------
     * 
     * 
     * The function print_region_probabilities(const RegionState& region) is not fundamentally necessary for the core functionality of the RegionGenerator class.
     * However, it can be a useful utility function for debugging and verification purposes.
     * It allows developers to easily inspect the base probabilities associated with a given region, ensuring that the probabilities align with the expected values based on the region's type and target content.
     * This can help in validating that the region generation logic is functioning correctly and that the probabilities are being calculated as intended.
     * 
     * 
     * --------------------------------------------------------------------------
     */

     void print_region_probabilities(const RegionState& region);
};

#endif // REGION_GENERATOR_H