#!/usr/bin/env python3
"""
build protein manifest from upstream segment_types and contacts_awsem data.
outputs: data/protein_manifest.csv with 37 Galpern proteins.
"""

import pandas as pd
from pathlib import Path

BASE = Path(__file__).resolve().parent.parent
DATA = BASE / "data"
UPSTREAM = DATA / "upstream"


def main():
    seg = pd.read_csv(UPSTREAM / "segment_types.csv")
    contacts = pd.read_csv(UPSTREAM / "contacts_awsem.csv")

    # filter to Galpern proteins only (protein_id starts with 'F')
    galpern_seg = seg[seg['protein_id'].str.startswith('F')].copy()

    rows = []
    for pid, grp in galpern_seg.groupby('protein_id'):
        pdb_id = grp['pdb_id'].iloc[0]
        n_res = grp['n_res'].iloc[0]
        n_foldons = len(grp)
        n_type_a = (grp['n_L2'] > 0).sum()
        n_type_b = (grp['n_L2'] == 0).sum()

        # count total L2 contacts from contacts_awsem
        pid_contacts = contacts[contacts['family_id'] == pid]
        n_l2 = (pid_contacts['layer'] == 'L2').sum()

        rows.append({
            'protein_id': pid,
            'pdb_id': pdb_id,
            'n_residues': n_res,
            'n_foldons': n_foldons,
            'n_type_a': n_type_a,
            'n_type_b': n_type_b,
            'n_l2_contacts': n_l2,
            'cds_status': 'pending',
        })

    manifest = pd.DataFrame(rows)
    # sort by protein_id naturally (F01, F02, ... F38)
    manifest['_sort'] = manifest['protein_id'].str.extract(r'(\d+)').astype(int)
    manifest = manifest.sort_values('_sort').drop(columns='_sort')

    out_path = DATA / "protein_manifest.csv"
    manifest.to_csv(out_path, index=False)
    print(f"wrote {len(manifest)} proteins to {out_path}")
    print(f"  type A foldons: {manifest['n_type_a'].sum()}")
    print(f"  type B foldons: {manifest['n_type_b'].sum()}")
    print(f"  total L2 contacts: {manifest['n_l2_contacts'].sum()}")


if __name__ == '__main__':
    main()
