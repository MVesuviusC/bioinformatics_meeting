import argparse
import os
import logging

def filter_vcf_by_chromosome(input_file, output_file, chromosome):
    """
    Filters a VCF file by a specific chromosome and writes the result to an output file.

    Args:
        input_file (str): Path to the input VCF file.
        output_file (str): Path to the output VCF file.
        chromosome (str): Chromosome to filter by.
    """
    try:
        with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
            for line in infile:
                if line.startswith("#"):
                    # Write header lines to the output file
                    outfile.write(line)
                else:
                    # Split the line to check the chromosome column
                    columns = line.split('\t')
                    if len(columns) > 0 and columns[0].lower() == chromosome.lower():
                        outfile.write(line)
    except FileNotFoundError:
        logging.error(f"Error: The file {input_file} does not exist.")
    except PermissionError:
        logging.error(f"Error: Permission denied for file {output_file}.")
    except Exception as e:
        logging.error(f"An unexpected error occurred: {e}")

def main():
    """
    Main function to parse arguments and filter the VCF file.
    """
    parser = argparse.ArgumentParser(description="Filter a VCF file by chromosome.")
    parser.add_argument("input_file", help="Path to the input VCF file.")
    parser.add_argument("output_file", help="Path to the output VCF file.")
    parser.add_argument("chromosome", help="Chromosome to filter by (e.g., 'chr1').")

    args = parser.parse_args()

    # Check if output file exists
    if os.path.exists(args.output_file):
        confirm = input(f"The file {args.output_file} already exists. Overwrite? (y/n): ")
        if confirm.lower() != 'y':
            print("Operation cancelled.")
            return

    logging.basicConfig(level=logging.INFO)
    logging.info("Starting the filtering process...")
    filter_vcf_by_chromosome(args.input_file, args.output_file, args.chromosome)
    logging.info("Filtering complete.")

if __name__ == "__main__":
    main()
