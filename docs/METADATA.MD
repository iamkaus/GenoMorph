enum class FeatureType {
    CDS,
    EXON,
    MRNA,
    GENE,
    INTRON,
    UTR5,
    UTR3,
    PROMOTER,
    ENHANCER,
    NCRNA,
    REPEAT_REGION,
    TELOMERE,
    CENTROMERE,
};

/**
 * --------------------------------------------------------------
 * NOTE: CODING REGION METADATA
 * --------------------------------------------------------------
 */

struct GeneIdentifiers {

private:
    std::string UNIQUE_ID;
    std::string GENE_NAME;
    std::string TRANSCRIPT_NAME;
};

struct GenomicCoordinates {

private:
    std::string     SCAFFOLD_ID;
    uint32_t        START_POSITION;
    uint32_t        END_POSITION;
};

struct Structure {

private:
    uint32_t TOTAL_EXONS;
    uint32_t EXON_LENGTH;
    uint32_t INTRON_LENGTH;
};

struct CodingSequence {

private:
    std::string OPENING_READING_FRAME_COORDINATES;
    std::string START_CODON;
    std::string END_CODON;
};

struct TranslationInfo {

private:
    std::string PROTEIN_PRODUCT_ID;
    std::string AMINO_ACID_SEQUENCE;
    std::string DOMAIN_ANNOTATION;
};

struct FunctionalAnnotation {
private:
    std::string     GENE_ONTOLOGY_TERMS;
    std::string     PATHWAY_ASSOCIATIONS;
    std::string     PROTEIN_PRODUCT_DESCRIPTION;
    uint32_t        ENZYME_COMMISSION_NUMBER;
    std::string     PROTEIN_DOMAINS;
};
/**
 * NOTE: GeneticVariation could be made into a struct needs an assurance whethere these are just constants or variable values
 */
enum class GeneticVariation {
    SYNONYMOUS_MUTATION,
    NON_SYNONYMOUS_MUTATION,
    FRAMESHIFTS,
    PREMATURE_STOP_CODONS
};

/**
 * --------------------------------------------------------
 * NOTE: NON-CODING REGION METADATA
 * --------------------------------------------------------
 */

enum class FunctionalNonCodingRNAS {
    MIRNA,
    TRNA,
    RRNA,
    SNORNA,
    INCRNA
};

enum class FunctionalRole {
    RNA_PROCESSING,
    GENE_REGULATION,
};

enum class PseudoGenes {
    INACTIVE,
    MUTATED_REMANANTS_OF_FORMER_GENES,
};

struct RegulatoryElements {

private:    
    struct ChromatinState {
        bool                IS_OPEN;
        bool                IS_CLOSED;

        /**
         * NOTE: HISTONE_MODIFICATION is highly variable process with multiple modification processes being possible under it
         */
        std::string         HISTONE_MODIFICATION;
    };

    struct TranscriptionFactorBindingSites {
        bool            BOUND_BY_PROTEINS;
        std::string     IDENTIFIED_BY_CHIP_SEQ;
    };

public:
    std::string                         PROMOTERS;
    std::string                         ENHANCERS;
    std::string                         SILENCERS;
    TranscriptionFactorBindingSites     TFBS;
    ChromatinState                      CSTATE;
};

struct NonCodingRNAGenes {

private:
    FunctionalNonCodingRNAS     TYPES;
    FunctionalRole              ROLE;
};

struct StructuralElements {
private:
    bool    IS_CENTROMERE;
    bool    IS_TELOMERE;
};

struct MobileGeneticElements {
private:
    std::string     TRANSPOSONS;
    std::string     ENDOGENOUS_RETROVIRUSES;
};

struct IntergenicRegions {

    /**
     * NOTE: IS_NON_FUNCTIONAL_DNA -> historically "junk" DNA
     * NOTE: CONSERVED_NON_CODING_ELEMENTS -> unknown function
     */
    bool        IS_NON_FUNCTIONAL_DNA;
    bool        CONSERVED_NON_CODING_ELEMENTS;
};

struct ConservationScores {

private:
    uint32_t    PHAST_CONS_SCORES;
    uint32_t    PHYLOP_SCORES;
};
