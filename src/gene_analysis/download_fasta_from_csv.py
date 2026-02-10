import csv
import argparse
import requests
import time
from Bio import SeqIO
from io import StringIO

EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

def fetch_fasta(accession, api_key=None):
    params = {
        "db": "nuccore",
        "id": accession,
        "rettype": "fasta",
        "retmode": "text",
    }
    if api_key:
        params["api_key"] = api_key

    r = requests.get(EFETCH_URL, params=params, timeout=30)
    r.raise_for_status()
    return r.text

def parse_fasta(text, accession):
    if ">" not in text:
        raise ValueError("No FASTA header returned")

    records = list(SeqIO.parse(StringIO(text), "fasta-blast"))
    if len(records) != 1:
        raise ValueError(f"Expected 1 FASTA record, got {len(records)}")

    return records[0]

def main(csv_file, output_fasta, api_key=None):
    records = []
    failed = []

    with open(csv_file, newline="") as f:
        reader = csv.DictReader(f)
        if "Accession" not in reader.fieldnames:
            raise ValueError("CSV must contain 'Accession' column")
        
        for row in reader:
            acc = row["Accession"].strip()
            print(f"Downloading {acc}...")
            

            try:
                fasta_text = fetch_fasta(acc, api_key)
                record = parse_fasta(fasta_text, acc)

                # Attach CSV metadata to header
                meta = [f"{k}={v}" for k, v in row.items() if k != "Accession"]
                if meta:
                    print(meta)
                    record.description += " | " + ";".join(meta)

                records.append(record)
                time.sleep(0.34)  # NCBI etiquette

            except Exception as e:
                print(f"  Skipping {acc}: {e}")
                failed.append(acc)

    if not records:
        raise RuntimeError("No sequences downloaded")

    SeqIO.write(records, output_fasta, "fasta")
    print(f"\nWrote {len(records)} sequences to {output_fasta}")

    if failed:
        print(f"⚠️  Failed accessions ({len(failed)}): {', '.join(failed)}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Download FASTA sequences from NCBI using accessions in CSV"
    )
    parser.add_argument("csv", help="Input CSV with Accession column")
    parser.add_argument("output", help="Output FASTA file")
    parser.add_argument("--api-key", help="NCBI API key (optional)")

    args = parser.parse_args()
    main(args.csv, args.output, args.api_key)
