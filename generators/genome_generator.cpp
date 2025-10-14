#include "genome_generator.h"
#include "regionGenerator.h"
#include <iostream>
#include <sstream>
#include <random>
#include <ctime>
#include <fstream>

/**
 * @brief Constructs a GenomeGenerator with optional base probabilities.
 *
 * Seeds the random number generator using the system clock.
 */
GenomeGenerator::GenomeGenerator()
{
    rng.seed(static_cast<unsigned>(time(0))); // Seed with system clock
}

/**
 * @brief Generates a single base according to weighted probabilities.
 * @param region The RegionState providing context for base generation.
 * @return A character representing the base ('A', 'T', 'C', or 'G').
 */
char GenomeGenerator::generate_base(RegionState region) {

    /**
     * @brief Determines base probabilities based on the region type.
     * Uses RegionGenerator to get region-specific probabilities.
     * Defaults to equal probabilities if region type is unknown.
     * This allows for biologically relevant base generation.
     */
    RegionGenerator regionGen;

    /**
     * @note DEBUG: The following print statement is for debugging purposes.
     */

    // std::cout << "Generating base for region type: " << region.region_type << " with target content: " << region.target_content << " and length: " << region.target_length << std::endl;

    /**
     * @note The regionBasedBaseProbabilities function returns an array of probabilities for A, T, C, and G.
     * These probabilities are then used to generate a base according to the specified weights.
     * The sum of these probabilities should ideally be 1.0.
     * This design allows for easy modification of base probabilities in the future.
     * Overall, this function encapsulates the logic for region-specific base probability determination.
     */
    auto probs = regionGen.regionBasedBaseProbabilities(region);
    regionGen.print_region_probabilities(region); // Debug: Print region probabilities

    double A_PROB = probs[0];
    double T_PROB = probs[1];
    double C_PROB = probs[2];
    double G_PROB = probs[3];

    std::uniform_real_distribution<double> dist(0.0, 1.0);
    double r = dist(rng);

    if ( r < A_PROB ) return 'A';
    else if ( r < A_PROB + T_PROB ) return 'T';
    else if ( r < A_PROB + T_PROB + G_PROB ) return 'G';
    else return 'C';
}

/**
 * @brief Generates a DNA sequence of given length.
 * @param length Number of bases in the sequence.
 * @param total_generated Length of genome generated so far.
 * @param filename Name of the output RTF file.
 * @return A vector of BaseInfo structs representing the sequence.
 * Each BaseInfo has default metadata values.
 */
std::vector<BaseInfo> GenomeGenerator::generate_sequence(size_t total_generated, size_t length, std::string &filename) {
    std::vector<BaseInfo> sequence;
    sequence.reserve(length);

    /**
     * @brief Uses RegionGenerator to create regions and generate bases accordingly.
     * Continues generating regions until the total length is reached.
     * Each base is generated based on the properties of the current region.
     * Metadata for each base is initialized to default values.
     * The function ensures that the generated sequence adheres to the specified length.
     * The use of RegionGenerator allows for realistic variation in the sequence structure.
     * The generated sequence is returned as a vector of BaseInfo structs.
     * This approach simulates the complexity of real genomic sequences.
    */
    RegionGenerator regionGen;

    /**
     * @brief Main loop to generate the sequence.
     * Continues until the total_generated length matches the desired length.
     * For each region, bases are generated and added to the sequence.
     * The current_content of the region is updated based on the bases generated.
     * Metadata such as coding_region, repair_efficiency, methylation, and chromatin_access are set for each base.
     * The loop ensures that the sequence is built incrementally, respecting region boundaries.
     * This method provides a structured way to create complex sequences with varying properties.
     * The final sequence is a combination of multiple regions, each with its own characteristics.
     * The function returns the complete sequence once the desired length is achieved.
     * This design allows for future enhancements, such as adding more metadata or region types.
     * Overall, this function encapsulates the logic for generating a biologically relevant DNA sequence.
    */
    while (total_generated < length) {
        RegionState region = regionGen.createRegion(total_generated, length);

        /**
         * @note DEBUG: The following print statement is for debugging purposes.
         */
        std::cout << "[Region Info] Type: " << region.region_type 
                  << " | Target GC Content: " << region.target_content 
                  << " | Target AT Content: " << (1.0 - region.target_content)
                  << " | Length: " << region.target_length << std::endl;

        for (size_t i = 0; i < region.target_length; ++i) {
            BaseInfo base;
            base.base = generate_base(region);

            /**
             * @brief Updates region content tracking.
             * Increments current_content if the base is G or C.
             * This helps monitor progress towards target_content.
            */
            if ( base.base == 'G' || base.base == 'c' ) {
                region.current_content += 1.0 / region.target_length;
            }

            /**
             * @brief Initializes base metadata.
             * Sets coding_region based on region type.
             * Other metadata are set to default values.
             * This metadata can be used for mutation simulations later.
             * The design allows for easy modification of metadata initialization in the future.
             * Overall, this ensures each base carries relevant biological information.
            */
            base.coding_region = (region.region_type == "coding");
            base.repair_efficiency = 1.0;
            base.methylation = 0.0;
            base.chromatin_access = 0.0;

            /**
             * @brief Adds the generated base to the sequence.
             * Increments the current_length_generated of the region.
             * This continues until the region's target_length is reached.
             * The loop then proceeds to generate the next region if needed.
             * The sequence is built incrementally, ensuring all bases are accounted for.
             * The final sequence reflects the combined properties of all generated regions.
            */
            sequence.push_back(base);
            region.current_length_generated++;
        }

        /**
         * @brief Updates the total length generated so far.
         * This is used to determine when the overall sequence generation is complete.
         * The loop continues until the desired length is achieved.
         * This ensures the function adheres to the specified length constraint.
         * The design allows for flexibility in region sizes and types.
         * Overall, this mechanism ensures a coherent and complete sequence generation process.
        */
        total_generated += region.target_length;
    }
    rtfGenome(sequence, filename); // Export the generated sequence to an RTF file
    return sequence;
}


/**
 * @brief Generates the complementary strand for a given vector of BaseInfo elements.
 * @param original The original single-stranded DNA sequence (5' → 3') as a vector of BaseInfo.
 * @return The complementary strand (3' → 5') as a vector of BaseInfo.
 */

 std::vector<BaseInfo> GenomeGenerator::ComplementaryStrand(const std::vector<BaseInfo> &original) {
    std::vector<BaseInfo> complementary_base;
    complementary_base.reserve(original.size());

    for (const auto& base_info: original) {
        BaseInfo comp_base = base_info;         // Copy metadata

        switch(base_info.base) {
            case 'A':
                comp_base.base = 'T';
                break;
            case 'T':
                comp_base.base = 'A';
                break;
            case 'C':
                comp_base.base = 'G';
                break;
            case 'G':
                comp_base.base = 'C';
                break;
            default:
                comp_base.base = 'N';          // Unknown base
                
            }
            /**
             * @note Metadata such as repair_efficiency, methylation, coding_region, and chromatin_access
             * are copied from the original base_info to the complementary base.
             * The properties of the complementary strand can be adjusted here if needed.
             * For example, methylation patterns might differ between strands in biological contexts.
             * However, for simplicity, we are keeping them the same in this implementation.
            */
        complementary_base.push_back(comp_base);
    }
    return complementary_base;
 }

/**
 * @brief Utility function to export the sequence to an RTF file.
 * @param sequence Vector of BaseInfo representing the sequence.
 * @param filename Name of the output RTF file.
 * Exports the sequence in a simple RTF format for easy viewing.
 */

 void GenomeGenerator::rtfGenome(const std::vector<BaseInfo> &sequence, const std::string &filename) {

    std::ofstream outfile(filename);
    if (!outfile.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return; 
    }

    outfile << "{\\rtf1\\ansi\\deff0\n";
    outfile << "{\\fonttbl{\\f0 Courier New;}}\n";

    outfile << "{\\colortbl;"
            << "\\red255\\green0\\blue0;"           // A → Red
            << "\\red0\\green0\\blue255;"           // T → Blue
            << "\\red0\\green200\\blue0;"           // G → Green
            << "\\red255\\green165\\blue0;}\n";     // C → Orange

    outfile << "\\f0\\fs20\n";                      // Set font and size
    outfile << "Colored DNA Sequence:\\line\n";

    const size_t wrap = 80;
    for (size_t i = 0; i < sequence.size(); ++i) {
        char base = sequence[i].base;

        switch (base) {
            case 'A': outfile << "\\cf1 " << base; break; // red
            case 'T': outfile << "\\cf2 " << base; break; // blue
            case 'G': outfile << "\\cf3 " << base; break; // green
            case 'C': outfile << "\\cf4 " << base; break; // orange
            default:  outfile << base; break;
        }

        if ((i + 1) % wrap == 0) outfile << "\\line\n";
    }

    outfile << "\\cf0\\line\n}\n";
    outfile.close();

    std::cout << "✅ Colored DNA sequence successfully written to " << filename << "\n";
 }

/**
 * @brief Prints the sequence to the console as a string of bases.
 * @param sequence Vector of BaseInfo representing the sequence.
 */
void GenomeGenerator::print_sequence(const std::vector<BaseInfo> &sequence) {
    for (const auto &b : sequence) {
        std::cout << b.base;
    }
    std::cout << std::endl;
}

/**
 * @brief Utility function to print both strands as a "double helix" representation.
 * @param strand1 The first strand (5' → 3').
 * @param strand2 The complementary strand (3' → 5').
 * Prints the strands in a paired format.
 */

 void GenomeGenerator::printDoubleHelix(const std::vector<BaseInfo> &strand1, const std::vector<BaseInfo> &strand2) {
     std::cout << "\n=== Double Helix Representation ===\n";
    for (size_t i = 0; i < strand1.size(); ++i) {
        std::cout << strand1[i].base << " - " << strand2[i].base << std::endl;
    }
 }

/**
 * @brief Converts the sequence to a string.
 * @param sequence Vector of BaseInfo representing the sequence.
 * @return std::string containing the sequence of bases.
 */

[[deprecated("Warning: This function is yet to be tested. Use with caution.")]]
std::string GenomeGenerator::sequence_to_string(const std::vector<BaseInfo> &sequence) {
    std::ostringstream oss;
    for (const auto &b : sequence) {
        oss << b.base;
    }
    return oss.str();
}
