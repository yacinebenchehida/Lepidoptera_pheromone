import re
import argparse

def filter_lines_by_species(input_file, output_file):
    # Regular expression to extract genus_species from a line
    species_pattern = re.compile(r'([A-Z][a-z]+_[a-z]+)')  # Matches "Genus_species"

    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            # Find all matches of genus_species in the current line
            species_found = species_pattern.findall(line)
            
            # Convert list to a set to get unique genus_species
            unique_species = set(species_found)
            
            # Check if there are at least two different genus_species
            if len(unique_species) >= 2:
                outfile.write(line)  # Keep the line in the output

def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Filter lines with at least two unique genus_species entries.")
    parser.add_argument("input_file", help="Path to the input file containing sequence data.")
    parser.add_argument("output_file", help="Path to the output file to save filtered results.")
    args = parser.parse_args()

    # Run the filtering function
    filter_lines_by_species(args.input_file, args.output_file)
    print(f"Filtered lines saved to: {args.output_file}")

if __name__ == "__main__":
    main()
