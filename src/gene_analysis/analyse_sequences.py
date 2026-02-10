import math
import re
import argparse
from datetime import datetime
from itertools import combinations
from io import StringIO

from Bio import AlignIO
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment

# ---------------------------------------------------------
# 1. Configuration & Constants
# ---------------------------------------------------------
DATE_RE = re.compile(r"Collection_Date=(\d{4}-\d{2}-\d{2})")

# HA1/HA2 Boundaries for H1N1/H3N2 (Approximate based on 1701bp)
HA1_START, HA1_END = 0, 987
HA2_START = 987 

# ---------------------------------------------------------
# 2. Core Metric Functions
# ---------------------------------------------------------

def calculate_shannon_entropy(alignment):
    """Calculates mean Shannon Entropy across all columns."""
    L = alignment.get_alignment_length()
    total_entropy = 0.0
    valid_cols = 0

    for i in range(L):
        col = [res for res in alignment[:, i] if res != "-"]
        if len(col) < 2: continue  # Skip low-coverage columns
        
        counts = {}
        for res in col:
            counts[res] = counts.get(res, 0) + 1
        
        entropy = 0.0
        for count in counts.values():
            p = count / len(col)
            entropy -= p * math.log2(p)
        
        total_entropy += entropy
        valid_cols += 1
        
    return total_entropy / valid_cols if valid_cols > 0 else 0.0

def percent_variable_sites(alignment):
    """Calculates the percentage of sites with more than one allele."""
    L = alignment.get_alignment_length()
    variable = 0
    valid_cols = 0

    for i in range(L):
        col = set([res for res in alignment[:, i] if res != "-"])
        if len(col) < 1: continue
        valid_cols += 1
        if len(col) > 1:
            variable += 1
            
    return (variable / valid_cols) * 100 if valid_cols > 0 else 0.0

def mean_identity(alignment):
    """Calculates average pairwise identity."""
    identities = []
    # If alignment is huge, we'd sample, but for ~500 seqs, we can do it
    seqs = [str(rec.seq) for rec in alignment]
    
    # Simple sampling if too many pairs to keep it fast
    import random
    pairs = list(combinations(range(len(seqs)), 2))
    if len(pairs) > 10000:
        pairs = random.sample(pairs, 10000)

    for i, j in pairs:
        s1, s2 = seqs[i], seqs[j]
        matches = sum(1 for a, b in zip(s1, s2) if a == b and a != "-" and b != "-")
        valid = sum(1 for a, b in zip(s1, s2) if a != "-" and b != "-")
        if valid > 0:
            identities.append(matches / valid)
            
    return (sum(identities) / len(identities)) * 100 if identities else 0.0

# ---------------------------------------------------------
# 3. Processing Functions
# ---------------------------------------------------------

def translate_alignment(nt_aln):
    """Translates a nucleotide alignment to Amino Acids codon-wise."""
    aa_recs = []
    for rec in nt_aln:
        aa_seq = []
        nt_str = str(rec.seq)
        for i in range(0, len(nt_str) - 2, 3):
            codon = nt_str[i:i+3]
            if "-" in codon:
                aa_seq.append("-")
            else:
                aa_seq.append(str(Seq(codon).translate()))
        
        new_rec = rec[:len(aa_seq)]
        new_rec.seq = Seq("".join(aa_seq))
        aa_recs.append(new_rec)
    return MultipleSeqAlignment(aa_recs)

def filter_by_date(alignment, start_str, end_str):
    """Filters alignment records by collection date."""
    start = datetime.fromisoformat(start_str)
    end = datetime.fromisoformat(end_str)
    
    kept = []
    for rec in alignment:
        m = DATE_RE.search(rec.description)
        if m:
            dt = datetime.fromisoformat(m.group(1))
            if start <= dt <= end:
                kept.append(rec)
    return MultipleSeqAlignment(kept) if len(kept) > 1 else None

# ---------------------------------------------------------
# 4. Main Script Execution
# ---------------------------------------------------------

def run_analysis(fasta_path, gene_name):
    print(f"\n{'='*40}")
    print(f"ANALYSIS FOR {gene_name}: {fasta_path}")
    print(f"{'='*40}")
    
    aln = AlignIO.read(fasta_path, "fasta")
    
    # Overall Stats
    print(f"Overall Identity: {mean_identity(aln):.2f}%")
    print(f"Overall Variable Sites: {percent_variable_sites(aln):.2f}%")
    print(f"Overall Mean Entropy: {calculate_shannon_entropy(aln):.4f}")

    # Temporal Windows
    windows = {
        "Early (Sept-Nov)": ("2024-09-01", "2024-11-30"),
        "Mid   (Nov-Jan)": ("2024-11-01", "2025-01-31"),
        "Late  (Jan-Mar)": ("2025-01-01", "2025-03-31")
    }

    for label, (s, e) in windows.items():
        sub_aln = filter_by_date(aln, s, e)
        if sub_aln:
            print(f"\n[{label}] (N={len(sub_aln)})")
            print(f"  Identity: {mean_identity(sub_aln):.2f}%")
            print(f"  Entropy:  {calculate_shannon_entropy(sub_aln):.4f}")
        else:
            print(f"\n[{label}]: No data.")

    # Domain Analysis for HA
    if "HA" in gene_name.upper():
        print(f"\n--- HA Domain Analysis (Full Season) ---")
        ha1 = aln[:, HA1_START:HA1_END]
        ha2 = aln[:, HA2_START:]
        
        print(f"HA1 (Head)  Entropy: {calculate_shannon_entropy(ha1):.4f}")
        print(f"HA2 (Stalk) Entropy: {calculate_shannon_entropy(ha2):.4f}")
        
        # AA Level
        print(f"\n--- Amino Acid Level Metrics ---")
        aa_aln = translate_alignment(aln)
        aa_ha1 = aa_aln[:, :HA1_END//3]
        aa_ha2 = aa_aln[:, HA1_END//3:]
        
        print(f"AA HA1 Variable Sites: {percent_variable_sites(aa_ha1):.2f}%")
        print(f"AA HA2 Variable Sites: {percent_variable_sites(aa_ha2):.2f}%")
        print(f"AA HA1 Entropy: {calculate_shannon_entropy(aa_ha1):.4f}")
        print(f"AA HA2 Entropy: {calculate_shannon_entropy(aa_ha2):.4f}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Full California Flu Season Analysis")
    parser.add_argument("--ha", help="Path to aligned curated HA file")
    parser.add_argument("--np", help="Path to aligned curated NP file")
    args = parser.parse_args()

    if args.ha: run_analysis(args.ha, "HA (Segment 4)")
    if args.np: run_analysis(args.np, "NP (Segment 5)")