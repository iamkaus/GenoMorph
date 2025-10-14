#include "regionGenerator.h"
#include "genome_generator.h"

/**
 * @brief Main function demonstrating genome generation and output.
 * Generates a DNA sequence, its complementary strand, and prints both.
 * Utilizes GenomeGenerator and RegionGenerator classes.
 * The generated sequence includes metadata for each base.
 * The complementary strand is derived from the original sequence.
 * Both strands are printed in a double helix format.
*/
int main() {

    /**
     * @brief Demonstrates the usage of GenomeGenerator.
     * Generates a sequence, its complementary strand, and prints both.
     * The sequence length is set to 1000 bases for this example.
     * The GenomeGenerator handles base probabilities and metadata initialization.
     * The output includes the original sequence, complementary strand, and a double helix representation.
    */
    GenomeGenerator generator;                                                            // A, T, C, G probabilities
    /**
     * @brief Specifies the output RTF file name.
     * The file will contain the generated DNA sequence in a colored format.
     * The filename can be modified as needed.
     * The RTF format allows for easy viewing of the sequence with color coding for different bases.
     * The file will be created in the "results" directory.
     * Ensure the directory exists before running the program to avoid file write errors.
     * The generated file can be opened with any RTF-compatible viewer.
     * This feature enhances the usability of the generated sequence data.
     * The filename variable is passed to the generate_sequence method.
     * This setup allows for easy customization of the output file location and name.
     * The chosen name "genome.rtf" reflects the content of the file.
     * The use of RTF format is a simple way to visualize the sequence without complex software.
     */
    std::string filename = "results/genome.rtf";                                       // Output RTF file name

    /**
     * @brief Generates a DNA sequence of specified length.
     * The length is set to 1000 bases in this example.
     */
    auto sequence = generator.generate_sequence(0, 10000, filename);                                // Generate sequence of length 50

    /**
     * @brief Generates the complementary strand.
     * Uses the ComplementaryStrand method of GenomeGenerator.
     */
    auto complementary = generator.ComplementaryStrand(sequence);                         // Generate complementary strand

    /**
     * @brief Converts the sequence to a string and prints it.
     * Uses the sequence_to_string method of GenomeGenerator.
     * This string representation can be useful for further analysis or output.
     * The conversion ensures that only the base characters are included, omitting metadata.
     * The resulting string is printed to the console.
     * This step demonstrates the utility of having a string representation of the sequence.
     * It can be easily integrated with other tools or libraries that operate on string data.
     * Overall, this enhances the usability of the generated sequence.
     * The string can also be saved to a file or used in downstream applications.
     * This final step completes the demonstration of the GenomeGenerator's capabilities.
     */

     /**
      * Note: The following line is commented out and not tested to avoid excessive console output.
      * Uncomment it if you wish to see the string representation of the sequence.
      * Be cautious as printing long sequences can clutter the console.
      */
    
    // std::string seq_str = generator.sequence_to_string(sequence);                         // Convert sequence to string
    // std::cout << "Sequence as string: " << seq_str << std::endl;                       // Print the string representation
    return 0;
}
