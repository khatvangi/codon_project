#!/usr/bin/env python3
"""
02_fetch_cds.py — fetch coding DNA sequences for all proteins in the manifest.

reads data/protein_manifest.csv, fetches CDS for each protein with
cds_status='pending', saves FASTA files to data/cds/, and updates manifest.
"""

import sys
import pandas as pd
from pathlib import Path

BASE = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(BASE / "src"))

from cds_mapping import fetch_cds_for_protein

MANIFEST = BASE / "data" / "protein_manifest.csv"
CDS_DIR = BASE / "data" / "cds"
CDS_DIR.mkdir(parents=True, exist_ok=True)


def main():
    df = pd.read_csv(MANIFEST)
    print(f"loaded manifest: {len(df)} proteins")

    # track results
    succeeded = []
    failed = []

    pending = df[df['cds_status'] == 'pending']
    print(f"  {len(pending)} with cds_status='pending'\n")

    for idx, row in pending.iterrows():
        protein_id = row['protein_id']
        pdb_id = row['pdb_id']
        print(f"[{protein_id}] fetching CDS for {pdb_id}...")

        try:
            cds_seq, prot_seq, organism = fetch_cds_for_protein(pdb_id)
        except Exception as e:
            print(f"  ERROR: {e}")
            df.at[idx, 'cds_status'] = 'failed'
            failed.append((protein_id, pdb_id, str(e)))
            continue

        if cds_seq is None:
            print(f"  FAILED: no CDS found")
            df.at[idx, 'cds_status'] = 'failed'
            failed.append((protein_id, pdb_id, "no CDS found"))
            continue

        # validate: CDS length should be multiple of 3
        if len(cds_seq) % 3 != 0:
            print(f"  WARNING: CDS length {len(cds_seq)} not divisible by 3, trimming")
            cds_seq = cds_seq[:len(cds_seq) - (len(cds_seq) % 3)]

        # save as FASTA
        fasta_path = CDS_DIR / f"{protein_id}.fasta"
        organism_clean = organism.replace('|', '_') if organism else 'unknown'
        header = f">{protein_id}|{pdb_id}|{organism_clean}"
        with open(fasta_path, 'w') as f:
            f.write(f"{header}\n{cds_seq}\n")

        n_codons = len(cds_seq) // 3
        print(f"  OK: {len(cds_seq)} nt ({n_codons} codons), organism={organism}")

        df.at[idx, 'cds_status'] = 'success'
        succeeded.append((protein_id, pdb_id, organism))

    # save updated manifest
    df.to_csv(MANIFEST, index=False)
    print(f"\n{'='*60}")
    print(f"SUMMARY: {len(succeeded)} succeeded, {len(failed)} failed")
    print(f"manifest saved to {MANIFEST}")

    if succeeded:
        print(f"\nsucceeded ({len(succeeded)}):")
        for pid, pdb, org in succeeded:
            print(f"  {pid} ({pdb}) — {org}")

    if failed:
        print(f"\nfailed ({len(failed)}):")
        for pid, pdb, reason in failed:
            print(f"  {pid} ({pdb}) — {reason}")


if __name__ == '__main__':
    main()
