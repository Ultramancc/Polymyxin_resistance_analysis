import os
from Bio import SeqIO
from Bio.Seq import Seq
import subprocess
import argparse
import sys

def translate_sequences(input_file, protein_file):
    """
    Translates nucleotide sequences to protein sequences.
    """
    with open(protein_file, "w") as output_handle:
        for record in SeqIO.parse(input_file, "fasta"):
            try:
                # Translate the sequence
                protein_seq = record.seq.translate(table=11)  # Standard Table 11 for bacteria/plastids
                # Write the protein sequence to the output file
                output_handle.write(f">{record.id}\n{str(protein_seq)}\n")
            except Exception as e:
                print(f"Error translating {record.id}: {e}")

def run_muscle(input_file, output_file):
    """
    Runs MUSCLE to generate a multiple sequence alignment using MUSCLE v5 syntax.
    """
    muscle_exe = "muscle"  # Ensure MUSCLE is installed and in your PATH
    muscle_cline = [muscle_exe, "-align", input_file, "-output", output_file]

    try:
        print(f"Running MUSCLE on {input_file}...")
        subprocess.run(muscle_cline, check=True)
        print(f"Alignment saved to {output_file}")
    except subprocess.CalledProcessError as e:
        print(f"Error running MUSCLE: {e}")
        raise
    except FileNotFoundError:
        print("Error: MUSCLE is not installed or not in your PATH.")
        raise

def find_mutations(alignment_file, mutation_table_file):
    """
    Compares aligned sequences to the reference sequence and identifies mutations.
    Outputs mutations in the format S238H.
    """
    # Parse the alignment
    alignment = list(SeqIO.parse(alignment_file, "fasta"))
    reference_seq = None

    # Identify the reference sequence
    for record in alignment:
        if record.id == "reference":
            reference_seq = str(record.seq)
            break

    # Exit with an error if no reference sequence is found
    if reference_seq is None:
        error_msg = f"Error: No reference sequence with ID 'reference' found in the alignment file '{alignment_file}'."
        print(error_msg)
        sys.exit(1)

    # Find mutations for each sample
    with open(mutation_table_file, "w") as output_handle:
        output_handle.write("Sample\tMutation\n")  # Header line
        for record in alignment:
            if record.id == "reference":
                continue  # Skip the reference sequence
            sample_id = record.id
            sample_seq = str(record.seq)

            mutations = []
            for i, (ref_residue, sample_residue) in enumerate(zip(reference_seq, sample_seq)):
                if ref_residue != sample_residue and ref_residue != "-" and sample_residue != "-":
                    mutation = f"{ref_residue}{i+1}{sample_residue}"  # Mutation format S238H
                    mutations.append(mutation)

            # Write mutations to the mutation table
            output_handle.write(f"{sample_id}\t{','.join(mutations)}\n")

    print(f"Mutation table saved to {mutation_table_file}")

if __name__ == "__main__":
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Translate nucleotide sequences to proteins, align them, and find mutations.")
    parser.add_argument("-i", "--input", required=True, help="Input nucleotide FASTA file (e.g., sequence.all.fna)")
    parser.add_argument("-o", "--output", required=True, help="Output alignment file (e.g., protein.aln.fna)")
    parser.add_argument("-t", "--table", required=True, help="Output mutation table file (e.g., mutations.tsv)")
    
    args = parser.parse_args()

    # Define intermediate protein file
    protein_file = "protein.all.fna"

    # Step 1: Translate nucleotide sequences to protein sequences
    print("Translating sequences...")
    translate_sequences(args.input, protein_file)

    # Step 2: Run MUSCLE to align protein sequences
    print("Generating alignment...")
    run_muscle(protein_file, args.output)

    # Step 3: Identify mutations
    print("Identifying mutations...")
    find_mutations(args.output, args.table)