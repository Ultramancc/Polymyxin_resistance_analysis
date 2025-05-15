import argparse


def parse_coverage(coverage):
    """
    Parses the coverage column (start-end/length) to extract start, end, and length.
    """
    range_part, length = coverage.split("/")
    start, end = map(int, range_part.split("-"))
    length = int(length)
    return start, end, length


def process_files(is_file, gene_file, output_file):
    """
    Processes the IS blastn results file and the gene output file to determine
    "before" and "after" insertion sites, based on the specified conditions.
    The results are written to the specified output file.
    """
    # Read the gene file into a list of dictionaries
    gene_data = []
    with open(gene_file, "r") as gf:
        header = gf.readline().strip().split("\t")  # Read and split the header
        for line in gf:
            values = line.strip().split("\t")
            gene_data.append(dict(zip(header, values)))  # Create a dictionary for each row

    # Open output file for writing
    with open(output_file, "w") as outf:
        # Write the header
        outf.write("File\tSEQUENCE\tInsertion_type\tIS\tInsertion_site\tGENE_name\n")

        # Process the IS results
        with open(is_file, "r") as isf:
            for line in isf:
                is_cols = line.strip().split("\t")
                is_sequence = is_cols[0]
                is_name = is_cols[1]  # IS name (column 2)
                is_start = int(is_cols[8])  # IS start position
                is_end = int(is_cols[9])    # IS end position

                for gene in gene_data:
                    if gene["SEQUENCE"] == is_sequence:
                        file_name = gene["#FILE"]
                        gene_start = int(gene["START"])
                        gene_end = int(gene["END"])
                        gene_name = gene["GENE"]
                        coverage = gene["COVERAGE"]

                        # Parse the coverage column
                        cov_start, cov_end, cov_length = parse_coverage(coverage)

                        # Determine "before" insertion site
                        before_insertion = None
                        if cov_start > 1:  # Check if coverage start > 1
                            threshold_before = gene_start - cov_start
                            if is_end > threshold_before:
                                before_insertion = "before"

                        # Determine "after" insertion site
                        after_insertion = None
                        if cov_end < cov_length:  # Check if coverage end < length
                            threshold_after = gene_end + (cov_length - cov_end)
                            if is_start < threshold_after:
                                after_insertion = "after"

                        # Write the results to the output file
                        if before_insertion:
                            outf.write(
                                f"{file_name}\t{is_sequence}\t{before_insertion}\t{is_name}\t{cov_start}\t{gene_name}\n"
                            )
                        if after_insertion:
                            outf.write(
                                f"{file_name}\t{is_sequence}\t{after_insertion}\t{is_name}\t{cov_end}\t{gene_name}\n"
                            )


if __name__ == "__main__":
    # Parse command-line arguments
    parser = argparse.ArgumentParser(
        description="Determine IS insertion sites relative to gene locations."
    )
    parser.add_argument(
        "-i", "--input", required=True, help="Input IS blastn results file (outfmt=6)."
    )
    parser.add_argument(
        "-g", "--gene", required=True, help="Input gene results file."
    )
    parser.add_argument(
        "-o", "--output", required=True, help="Output file to save the results."
    )

    args = parser.parse_args()

    # Process the files
    process_files(args.input, args.gene, args.output)
    print(f"Results saved to {args.output}")