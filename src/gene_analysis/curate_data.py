import argparse
from Bio import SeqIO
import os

def curate_influenza_data(input_file, output_file, target_len):
    kept = []
    atg_count = 0
    
    if not os.path.exists(input_file):
        print(f"Error: {input_file} not found.")
        return

    for record in SeqIO.parse(input_file, "fasta"):
        # Check for exact length match
        if len(record.seq) == target_len:
            kept.append(record)
            # Check if it starts with the start codon
            if str(record.seq).upper().startswith("ATG"):
                atg_count += 1
    
    if kept:
        SeqIO.write(kept, output_file, "fasta")
        print(f"--- Curation Summary for {input_file} ---")
        print(f"Target Length: {target_len} bp")
        print(f"Sequences Kept: {len(kept)} / {atg_count} start with ATG")
        print(f"Saved to: {output_file}\n")
    else:
        print(f" No sequences of length {target_len} found in {input_file}.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Filter FASTA sequences by exact length for high-quality genomic analysis."
    )
    parser.add_argument("input", help="Path to the raw input FASTA file")
    parser.add_argument("output", help="Path to save the curated output FASTA file")
    parser.add_argument("length", type=int, help="The target nucleotide length (e.g., 1701 for HA)")

    args = parser.parse_args()
    
    curate_influenza_data(args.input, args.output, args.length)