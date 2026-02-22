#pragma once

#include <string>
#include <random>
#include <array>
#include <optional>
#include <vector>

/**
 * @struct CodingMetaData
 * @brief struct with reading_frame as its member it described one of the three possible way to parse a nucleotide as a codon.
 * 
 * Codon -> a codon is a consecutive non-overlapping triplet of nucleotide.
 */

struct CodingMetaData { 
    /**
     * NODE: READING_FREAME -> can either be (1,2,3), (-1,-2,-3) or ORF (has a start codon [ATG] and an end codon [TAA,TAG, TGA])
     */
    int8_t reading_frame;
};

/**
 * @struct RegulatoryMetaData
 * @brief the struct refers to non-coding information embedded within genome that dictates [when, where and how much] a gene is "expressed". 
 */

struct RegulatoryMetaData { 
    /**
     * NOTE: ACCESSIBILITY -> between 0.0 - 1.0
     */
    double accessibility;
};

/**
 * @enum FeatureType
 * @brief enum FeatureType consists of integral constants to determine feature type of nucleotide bases.
 */

enum class FeatureType {

    coding,
    non_coding,
    regulatory,
    repeat
};

/**
 * @enum StrandInfo
 * @brief it provides information regarding a convention to specify which of the two strand in a double stranded DNA molecule contains a specific gene of sequence feature.
 */

enum class StrandInfo {plus, minus};

/**
 * NOTE: RegionPlan -> structural roadmap to DNA regions with common shared metadata
 */

struct RegionPlan {
    size_t       region_start_index;
    size_t       region_end_index;
    StrandInfo   strand;
    
    size_t RegionLength () const { return region_end_index - region_start_index + 1; }
};

/**
 * NOTE: BaseRegionInfo -> common shared metadata for base pairs
 */

struct BaseRegionInfo {
    FeatureType     type;
    RegionPlan      region_plan;
    double          AT_CONTENT =         0.0;
    double          GC_CONTENT =         0.0;
};

/**
 * @struct RegionInfo
 * @brief conveys information regarding region bases and whethere the region is coding or non-coding depending upon optional members.
 */

struct RegionInfo {
    BaseRegionInfo                                  base;
    std::optional<CodingMetaData>                   coding;
    std::optional<RegulatoryMetaData>               regulatory_meta_data;
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
     * @return A RegionInfo struct representing the newly created region.
     * Throws std::invalid_argument if genomeLength is less than 100.
    */
    RegionInfo createRegion(size_t currentGenomeLength, size_t genomeLength);

    /**
     * @brief Provides base probabilities based on the region type.
     * @param region The RegionState providing context for base probability determination.
     * @return An array of doubles representing the probabilities for A, T, C, and G respectively.
     * This function adjusts base probabilities to reflect the characteristics of the region.
     * For coding regions, higher GC content is favored, while non-coding regions favor AT content.
     * The returned probabilities can be used in base generation to ensure the sequence adheres to the region's properties.
     */
    std::array<double, 4> regionBasedBaseProbabilities(const RegionInfo &region);
};