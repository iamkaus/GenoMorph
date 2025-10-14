#ifndef GENOME_GENERATOR_H
#define GENOME_GENERATOR_H

#include "regionGenerator.h"

#include <vector>
#include <string>
#include <random>
#include <ctime>

/**
 * @struct BaseInfo
 * @brief Represents a single nucleotide base in a DNA sequence along with metadata for mutation simulation.
 *
 * This struct stores the identity of the base (A, T, C, G) and optional biological or experimental
 * attributes that can influence mutation probability or other sequence analyses.
 */
struct BaseInfo {
    /**
     * @brief The nucleotide base.
     * Possible values: 'A', 'T', 'C', 'G'.
     */
    char base;

    /**
     * @brief DNA repair efficiency at this base.
     * Value between 0.0 and 1.0 (default 1.0). Higher values mean more efficient repair,
     * reducing mutation likelihood.
     */
    double repair_efficiency;

    /**
     * @brief Methylation status at this base.
     * Value between 0.0 (unmethylated) and 1.0 (fully methylated). Default is 0.0.
     * Methylation can influence mutation probability, e.g., C→T transitions at CpG sites.
     * Methylation is a biochemical process that modifies DNA by adding a methyl group (–CH₃) to certain DNA bases, usually cytosine. It’s an important epigenetic mechanism, meaning it can regulate gene activity without changing the underlying DNA sequence.
     */
    double methylation;

    /**
     * @brief Whether this base is part of a coding region.
     * A coding region in DNA is the part of a gene that contains the instructions to make a protein. It’s essentially the “blueprint” for building the amino acid sequence of a protein.
     * Default is false. Useful for applying different mutation rules for coding vs non-coding DNA.
     */
    bool coding_region;

    /**
     * @brief Chromatin accessibility of this base.
     * Chromatin accessibility (or chromatin access) refers to how “open” or “closed” the DNA is inside the cell nucleus, which affects whether genes can be read (transcribed) into RNA.
     * Value between 0.0 (closed chromatin) and 1.0 (open chromatin). Default is 0.0.
     * Open chromatin may have higher exposure to mutagens.
     */
    double chromatin_access;
};

/**
 * @class GenomeGenerator
 * @brief Generates DNA sequences with optional weighted base probabilities and metadata.
 *
 * The GenomeGenerator class provides functionality to generate DNA sequences of arbitrary length,
 * store additional base-level metadata for mutation simulation, and export or print sequences.
 * Base probabilities can be weighted to simulate realistic GC content or other nucleotide distributions.
 */
class GenomeGenerator {
private:
    std::mt19937 rng; /**< Random number generator seeded with system clock. */

    /**
     * @brief Generates a single base according to weighted probabilities.
     * @param region The RegionState providing context for base generation.
     * @return A character representing the base ('A', 'T', 'C', or 'G').
     */
    char generate_base(RegionState region);

public:
    /**
     * @brief Constructs a GenomeGenerator with optional base probabilities.
     * @param A Probability of generating 'A' (default 0.25).
     * @param T Probability of generating 'T' (default 0.25).
     * @param C Probability of generating 'C' (default 0.25).
     * 
     * @param G Probability of generating 'G' (default 0.25).
     *
     * The sum of A, T, C, and G should ideally be 1.0. Probabilities are used for weighted random generation.
     */
    GenomeGenerator();

    /**
     * @brief Generates a DNA sequence of given length.
     * @param length The number of bases in the sequence.
     * @param total_generated Length of genome generated so far.
     * @param filename Name of the output RTF file.
     * @return A vector of BaseInfo structs representing the generated sequence.
     * Each BaseInfo contains the base and default metadata values (repair efficiency, methylation, coding region, chromatin accessibility).
     */
    std::vector<BaseInfo> generate_sequence(size_t currentGenomeLength, size_t length, std::string &filename);
    
    /**
     * @brief Generates the complementary strand for a given DNA sequence.
     * @param original The original single-stranded DNA sequence (5' → 3').
     * @return The complementary strand (3' → 5').
     */
    std::vector<BaseInfo> ComplementaryStrand(const std::vector<BaseInfo> &original); 

    /**
     * @brief Utility function to export the sequence to an RTF file.
     * @param sequence The vector of BaseInfo representing the sequence.
     * @param filename The name of the output RTF file.
     * Exports the sequence in a simple RTF format for easy viewing.
     */
    void rtfGenome(const std::vector<BaseInfo> &sequence, const std::string &filename);

    /**
     * @brief Utility function to print both strands as a "double helix" representation.
     * @param strand1 The first strand (5' → 3').
     * @param strand2 The complementary strand (3' → 5').
     * Prints the strands in a paired format.
     */
    void printDoubleHelix(const std::vector<BaseInfo> &strand1, const std::vector<BaseInfo> &strand2);


    /**
     * @brief Prints the sequence to the console as a string of bases.
     * @param sequence The vector of BaseInfo representing the sequence.
     */
    void print_sequence(const std::vector<BaseInfo> &sequence);

    /**
     * @brief Converts the sequence to a string.
     * @param sequence The vector of BaseInfo representing the sequence.
     * @return A std::string containing the sequence of bases.
     */

    [[deprecated("Warning: This function is yet to be tested. Use with cautionBe cautious as printing long sequences can clutter the console..")]]
    std::string sequence_to_string(const std::vector<BaseInfo> &sequence);
};


#endif // GENOME_GENERATOR_H
