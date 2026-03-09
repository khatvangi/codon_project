"""
cds_mapping.py — fetch coding DNA sequences for PDB proteins.

pipeline: PDB → SIFTS → UniProt → EMBL cross-ref → NCBI GenBank CDS.
ported from foldon_project/codon_layer_analysis.py with dual-cache support:
  - upstream caches (symlinked from foldon_project, read-only)
  - local data/cds/ cache (for new downloads)
"""

import json
import time
import requests
from pathlib import Path
from Bio import Entrez, SeqIO

BASE = Path(__file__).resolve().parent.parent

# upstream caches (read-only, from foldon_project via symlinks)
UPSTREAM_SIFTS = BASE / "data" / "upstream" / "uniprot_cache"
UPSTREAM_CDS = BASE / "data" / "upstream" / "cds_cache"

# local cache (new downloads go here)
LOCAL_CDS = BASE / "data" / "cds"
LOCAL_CDS.mkdir(parents=True, exist_ok=True)

Entrez.email = "research@cotrans-layer.org"
API_DELAY = 0.4


# ──────────────────────────────────────────────────────────────────────
# SIFTS: PDB → UniProt mapping
# ──────────────────────────────────────────────────────────────────────

def fetch_sifts_mapping(pdb_id):
    """
    fetch PDB → UniProt mapping via SIFTS API.
    checks upstream cache first, then local cache, then fetches from API.
    returns dict with keys: accession, chain, residue_map, n_mapped
    """
    # check upstream cache
    cache_file = UPSTREAM_SIFTS / f"{pdb_id}_sifts.json"
    if cache_file.exists():
        with open(cache_file) as f:
            return json.load(f)

    # check local cache
    local_file = LOCAL_CDS / f"{pdb_id}_sifts.json"
    if local_file.exists():
        with open(local_file) as f:
            return json.load(f)

    # fetch from PDBe API
    url = f"https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{pdb_id}"
    try:
        resp = requests.get(url, timeout=30)
        resp.raise_for_status()
        data = resp.json()
    except Exception as e:
        print(f"    SIFTS error for {pdb_id}: {e}")
        return None

    pdb_data = data.get(pdb_id, {}).get('UniProt', {})
    if not pdb_data:
        return None

    accession = list(pdb_data.keys())[0]
    mappings = pdb_data[accession].get('mappings', [])

    # build residue map from first chain encountered
    residue_map = {}
    chain_used = None
    for mapping in mappings:
        chain = mapping.get('chain_id', 'A')
        if chain_used is None:
            chain_used = chain
        if chain != chain_used:
            continue
        pdb_start = mapping.get('start', {}).get('residue_number')
        pdb_end = mapping.get('end', {}).get('residue_number')
        unp_start = mapping.get('unp_start')
        unp_end = mapping.get('unp_end')
        if all(v is not None for v in [pdb_start, pdb_end, unp_start, unp_end]):
            for offset in range(pdb_end - pdb_start + 1):
                residue_map[pdb_start + offset] = unp_start + offset

    result = {
        'accession': accession,
        'chain': chain_used,
        'residue_map': {str(k): v for k, v in residue_map.items()},
        'n_mapped': len(residue_map),
    }

    # save to local cache
    with open(local_file, 'w') as f:
        json.dump(result, f, indent=2)

    time.sleep(API_DELAY)
    return result


# ──────────────────────────────────────────────────────────────────────
# UniProt → EMBL cross-references
# ──────────────────────────────────────────────────────────────────────

def fetch_uniprot_cds_refs(accession):
    """
    fetch EMBL/GenBank CDS cross-references from UniProt.
    checks upstream cache, then local cache, then fetches from API.
    returns dict with keys: accession, organism, embl_refs
    """
    # check upstream cache
    cache_file = UPSTREAM_CDS / f"{accession}_cds_refs.json"
    if cache_file.exists():
        with open(cache_file) as f:
            return json.load(f)

    # check local cache
    local_file = LOCAL_CDS / f"{accession}_cds_refs.json"
    if local_file.exists():
        with open(local_file) as f:
            return json.load(f)

    # fetch from UniProt API
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.json"
    try:
        resp = requests.get(url, timeout=30)
        resp.raise_for_status()
        data = resp.json()
    except Exception as e:
        print(f"    UniProt API error for {accession}: {e}")
        return None

    # extract EMBL cross-references with CDS protein IDs
    embl_refs = []
    for xref in data.get('uniProtKBCrossReferences', []):
        if xref.get('database') == 'EMBL':
            nuc_id = xref.get('id')
            protein_id = None
            mol_type = None
            for prop in xref.get('properties', []):
                if prop.get('key') == 'ProteinId':
                    protein_id = prop.get('value')
                elif prop.get('key') == 'MoleculeType':
                    mol_type = prop.get('value')

            if nuc_id and protein_id and protein_id != '-':
                embl_refs.append({
                    'nucleotide_id': nuc_id,
                    'protein_id': protein_id,
                    'molecule_type': mol_type,
                })

    # grab organism name
    organism = data.get('organism', {}).get('scientificName', 'unknown')

    result = {
        'accession': accession,
        'organism': organism,
        'embl_refs': embl_refs,
    }

    # save to local cache
    with open(local_file, 'w') as f:
        json.dump(result, f, indent=2)

    time.sleep(API_DELAY)
    return result


# ──────────────────────────────────────────────────────────────────────
# NCBI GenBank CDS fetch
# ──────────────────────────────────────────────────────────────────────

def fetch_cds_from_ncbi(nucleotide_id, protein_id=None):
    """
    fetch CDS nucleotide sequence from NCBI GenBank.
    tries to find CDS feature matching protein_id, falls back to first CDS.
    checks upstream cache, then local cache, then fetches from API.
    returns dict with keys: nucleotide_id, cds_sequence, protein_sequence,
                            length_nt, length_aa, matched_protein_id
    """
    # check upstream cache
    cache_file = UPSTREAM_CDS / f"{nucleotide_id}_cds.json"
    if cache_file.exists():
        with open(cache_file) as f:
            return json.load(f)

    # check local cache
    local_file = LOCAL_CDS / f"{nucleotide_id}_cds.json"
    if local_file.exists():
        with open(local_file) as f:
            return json.load(f)

    # fetch from NCBI
    try:
        handle = Entrez.efetch(db="nucleotide", id=nucleotide_id,
                               rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        handle.close()
    except Exception as e:
        print(f"    NCBI error for {nucleotide_id}: {e}")
        return None

    # find CDS feature
    cds_seq = None
    protein_seq = None
    matched_protein_id = False

    for feature in record.features:
        if feature.type != 'CDS':
            continue
        qualifiers = feature.qualifiers
        prot_id = qualifiers.get('protein_id', [''])[0]

        try:
            cds_nt = str(feature.extract(record.seq))
            cds_prot = qualifiers.get('translation', [''])[0]
        except Exception:
            continue

        if not cds_nt or not cds_prot:
            continue

        # prefer CDS matching our protein_id
        if protein_id and protein_id.split('.')[0] in prot_id:
            cds_seq = cds_nt
            protein_seq = cds_prot
            matched_protein_id = True
            break
        elif cds_seq is None:
            cds_seq = cds_nt
            protein_seq = cds_prot

    # fallback: if entire record is an mRNA, use whole sequence
    if cds_seq is None and len(record.seq) % 3 == 0 and len(record.seq) > 60:
        from Bio.Seq import Seq
        try:
            test_prot = str(Seq(str(record.seq)).translate())
            if '*' not in test_prot[:-1]:  # no internal stops
                cds_seq = str(record.seq)
                protein_seq = test_prot.rstrip('*')
        except Exception:
            pass

    if cds_seq is None:
        return None

    result = {
        'nucleotide_id': nucleotide_id,
        'cds_sequence': cds_seq,
        'protein_sequence': protein_seq,
        'length_nt': len(cds_seq),
        'length_aa': len(protein_seq) if protein_seq else 0,
        'matched_protein_id': matched_protein_id,
    }

    # save to local cache
    with open(local_file, 'w') as f:
        json.dump(result, f, indent=2)

    time.sleep(API_DELAY)
    return result


# ──────────────────────────────────────────────────────────────────────
# full pipeline: PDB → CDS
# ──────────────────────────────────────────────────────────────────────

def fetch_cds_for_protein(pdb_id):
    """
    full pipeline: PDB → SIFTS → UniProt → EMBL → NCBI CDS.
    returns (cds_sequence, protein_sequence, organism) or (None, None, None).
    """
    # step 1: PDB → UniProt
    sifts = fetch_sifts_mapping(pdb_id)
    if not sifts:
        return None, None, None

    accession = sifts['accession']

    # step 2: UniProt → EMBL cross-refs
    refs = fetch_uniprot_cds_refs(accession)
    if not refs or not refs.get('embl_refs'):
        print(f"    no EMBL CDS refs for {accession}")
        return None, None, None

    organism = refs.get('organism', 'unknown')

    # step 3: try each EMBL ref until we get a CDS
    # prefer mRNA entries over genomic
    sorted_refs = sorted(refs['embl_refs'],
                         key=lambda x: (x.get('molecule_type', '') != 'mRNA',
                                        x['nucleotide_id']))

    for ref in sorted_refs[:5]:  # try up to 5
        cds = fetch_cds_from_ncbi(ref['nucleotide_id'], ref.get('protein_id'))
        if cds and cds.get('cds_sequence') and cds.get('protein_sequence'):
            return cds['cds_sequence'], cds['protein_sequence'], organism

    print(f"    no CDS found in {len(refs['embl_refs'])} EMBL refs")
    return None, None, None


# ──────────────────────────────────────────────────────────────────────
# sequence alignment: CDS protein → PDB residues
# ──────────────────────────────────────────────────────────────────────

def align_cds_to_pdb(cds_protein, pdb_sequence):
    """
    align CDS protein translation to PDB sequence.
    returns (mapping, n_matched):
      mapping: dict of pdb_res_index (0-based) → cds_res_index (0-based)
      n_matched: number of matched residues
    uses substring search first, then sliding window.
    """
    # exact substring: PDB is subsequence of CDS protein
    idx = cds_protein.find(pdb_sequence)
    if idx >= 0:
        return {i: idx + i for i in range(len(pdb_sequence))}, len(pdb_sequence)

    # exact substring: CDS protein is subsequence of PDB (rare)
    idx = pdb_sequence.find(cds_protein)
    if idx >= 0:
        return {idx + i: i for i in range(len(cds_protein))}, len(cds_protein)

    # fuzzy: sliding window with best match
    best_score = -1
    best_offset = 0
    pdb_len = len(pdb_sequence)
    cds_len = len(cds_protein)

    for offset in range(max(1, cds_len - pdb_len + 1)):
        match_len = min(pdb_len, cds_len - offset)
        score = sum(1 for i in range(match_len)
                    if pdb_sequence[i] == cds_protein[offset + i])
        if score > best_score:
            best_score = score
            best_offset = offset

    match_len = min(pdb_len, cds_len - best_offset)
    if match_len > 0 and best_score / match_len < 0.7:
        return None, 0  # too poor a match

    mapping = {}
    matched = 0
    for i in range(match_len):
        if pdb_sequence[i] == cds_protein[best_offset + i]:
            mapping[i] = best_offset + i
            matched += 1

    return mapping, matched


# ──────────────────────────────────────────────────────────────────────
# codon extraction
# ──────────────────────────────────────────────────────────────────────

def extract_codons(cds_sequence):
    """
    split a DNA string into a list of codon triplets (3-letter strings).
    any trailing nucleotides that don't form a complete triplet are dropped.
    """
    cds_sequence = cds_sequence.upper().replace('U', 'T')
    n = len(cds_sequence)
    return [cds_sequence[i:i+3] for i in range(0, n - (n % 3), 3)]
