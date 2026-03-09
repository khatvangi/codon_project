# Co-Translational Folding × Layer Architecture — Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Test whether L2 (intra-foldon) contacts form co-translationally before L3 (inter-foldon) contacts become geometrically possible, and whether this ordering is architecture-determined (invariant across translation rate regimes).

**Architecture:** O'Brien kinetic model applied per-foldon, with per-foldon contact order (CO) determining kf via Plaxco relation. Three rate regimes (uniform, tAI, shuffled-synonymous) × 3 tunnel lengths × 4 ku/kf ratios = 36 parameter combinations. Statistical tests across 37 proteins with 165 type A foldons.

**Tech Stack:** Python 3.11+, numpy, scipy, pandas, matplotlib, seaborn, requests, biopython

---

## Task 1: Create Directory Structure and Utils

**Files:**
- Create: `cotrans-layer/src/__init__.py`
- Create: `cotrans-layer/src/utils.py`
- Create: `cotrans-layer/tests/__init__.py`

**Step 1: Create directory structure**

```bash
cd /storage/kiran-stuff/codon_project
mkdir -p cotrans-layer/{data/{cds,upstream},src,scripts,results/{folding_curves,sensitivity,figures,summary},tests}
```

**Step 2: Create symlinks to upstream data**

```bash
cd /storage/kiran-stuff/codon_project/cotrans-layer/data/upstream
ln -s /storage/kiran-stuff/foldon_project/layer_results/segment_types.csv .
ln -s /storage/kiran-stuff/foldon_project/layer_results/contacts_awsem.csv .
ln -s /storage/kiran-stuff/foldon_project/paper_figures/tables/table_s1_all_proteins.csv .
ln -s /storage/kiran-stuff/foldon_project/uniprot_cache uniprot_cache
ln -s /storage/kiran-stuff/foldon_project/codon_results/cds_cache cds_cache
```

**Step 3: Write src/utils.py**

Contains:
- `CODON_TABLE`: standard genetic code dict (64 codons → 20 AA + stop)
- `AA_CODONS`: dict mapping each AA → list of synonymous codons
- `load_tai(path)`: load ecoli_tai_ws.tsv, return dict codon → tAI value
- `load_mu(path)`: load codon_error_rates.tsv, return dict codon → mu value
- `load_segment_types(path)`: load segment_types.csv, return DataFrame
- `load_contacts(path, family_id=None)`: load contacts_awsem.csv, optionally filter by family_id, return DataFrame

Data paths:
- tAI: `/storage/kiran-stuff/codon_project/ecoli_tai_ws.tsv`
- mu: `/storage/kiran-stuff/codon_project/codon_error_rates.tsv`
- segment_types: `data/upstream/segment_types.csv`
- contacts: `data/upstream/contacts_awsem.csv`

**Step 4: Write src/__init__.py and tests/__init__.py**

Empty files.

**Step 5: Commit**

```bash
git add cotrans-layer/
git commit -m "scaffold: cotrans-layer directory structure and utils"
```

---

## Task 2: Build Protein Manifest

**Files:**
- Create: `cotrans-layer/scripts/01_build_manifest.py`
- Create: `cotrans-layer/data/protein_manifest.csv`
- Test: manual inspection

**Step 1: Write 01_build_manifest.py**

Reads segment_types.csv and contacts_awsem.csv.
For each of the 37 Galpern proteins (protein_id starts with 'F'):
- Extract pdb_id, n_residues, n_foldons (= n_segments)
- Count type A segments (n_L2 > 0), type B segments (n_L2 = 0)
- Count total L2 contacts from contacts_awsem.csv

Output columns:
```
protein_id,pdb_id,n_residues,n_foldons,n_type_a,n_type_b,n_l2_contacts,cds_status
```

`cds_status` starts as 'pending' for all. Updated by script 02.

**Step 2: Run it**

```bash
cd /storage/kiran-stuff/codon_project/cotrans-layer
python scripts/01_build_manifest.py
```

Expected: 37 rows in `data/protein_manifest.csv`

**Step 3: Verify**

```bash
wc -l data/protein_manifest.csv  # should be 38 (37 + header)
head -5 data/protein_manifest.csv
```

**Step 4: Commit**

```bash
git add scripts/01_build_manifest.py data/protein_manifest.csv
git commit -m "feat: build protein manifest for 37 Galpern proteins"
```

---

## Task 3: Kinetic Model Core — Contact Order and P_folded

**Files:**
- Create: `cotrans-layer/src/kinetic_model.py`
- Create: `cotrans-layer/tests/test_kinetic_model.py`

**Step 1: Write failing tests in test_kinetic_model.py**

```python
import numpy as np
import pytest
from src.kinetic_model import compute_contact_order, compute_kf, cotrans_folding_probability

def test_contact_order_simple():
    """two contacts at known separations in a 10-residue foldon"""
    contacts = [(0, 5), (2, 8)]  # seq_sep = 5, 6
    co = compute_contact_order(contacts, foldon_length=10)
    expected = (1.0 / (10 * 2)) * (5 + 6)  # = 0.55
    assert abs(co - expected) < 1e-10

def test_contact_order_empty():
    """zero contacts → CO undefined, return None"""
    co = compute_contact_order([], foldon_length=10)
    assert co is None

def test_kf_from_plaxco():
    """known CO → known kf via Plaxco relation"""
    co = 0.2
    kf = compute_kf(co)
    # ln(kf) = -17.0 * 0.2 + 8.0 = 4.6
    expected = np.exp(4.6)
    assert abs(kf - expected) < 1e-6

def test_pfold_before_emergence():
    """chain length before emergence → P_folded = 0"""
    codon_times = np.full(100, 0.06)  # 100 codons, uniform
    P = cotrans_folding_probability(
        seg_start=10, seg_end=30, kf=100.0, ku=0.1,
        codon_times=codon_times, tunnel_length=30
    )
    # emergence at 30 + 30 = 60
    assert np.all(P[:61] == 0.0)

def test_pfold_approaches_equilibrium():
    """long time after emergence → P_folded → kf/(kf+ku)"""
    codon_times = np.full(200, 0.06)
    kf, ku = 100.0, 0.1
    P = cotrans_folding_probability(
        seg_start=10, seg_end=30, kf=kf, ku=ku,
        codon_times=codon_times, tunnel_length=30
    )
    P_eq = kf / (kf + ku)
    assert abs(P[-1] - P_eq) < 0.01  # within 1% of equilibrium

def test_pfold_monotonically_increasing():
    """P_folded should never decrease with chain length"""
    codon_times = np.full(150, 0.06)
    P = cotrans_folding_probability(
        seg_start=5, seg_end=20, kf=50.0, ku=0.05,
        codon_times=codon_times, tunnel_length=30
    )
    diffs = np.diff(P)
    assert np.all(diffs >= -1e-15)  # allow tiny floating-point noise
```

**Step 2: Run tests to verify they fail**

```bash
cd /storage/kiran-stuff/codon_project/cotrans-layer
python -m pytest tests/test_kinetic_model.py -v
```

Expected: all FAIL with ImportError

**Step 3: Write src/kinetic_model.py**

Functions:
- `compute_contact_order(contacts, foldon_length)`:
  - contacts: list of (i, j) residue pairs (0-indexed within the foldon)
  - Returns CO = (1 / (L * N)) * sum(|i-j|), or None if no contacts

- `compute_kf(contact_order)`:
  - Returns exp(-17.0 * CO + 8.0) (Plaxco relation, units: s^-1)

- `cotrans_folding_probability(seg_start, seg_end, kf, ku, codon_times, tunnel_length=30)`:
  - Returns array P_folded of length len(codon_times)
  - P_folded[L] = 0 for L < seg_end + tunnel_length
  - P_folded[L] = kf/(kf+ku) * (1 - exp(-(kf+ku) * t_avail)) for L >= emergence
  - t_avail = sum(codon_times[emergence:L+1])

**Step 4: Run tests to verify they pass**

```bash
python -m pytest tests/test_kinetic_model.py -v
```

Expected: all PASS

**Step 5: Commit**

```bash
git add src/kinetic_model.py tests/test_kinetic_model.py
git commit -m "feat: kinetic model core — contact order, kf, P_folded"
```

---

## Task 4: Contact Analysis — Per-Foldon CO from Real Data

**Files:**
- Create: `cotrans-layer/src/contact_analysis.py`
- Create: `cotrans-layer/tests/test_contact_analysis.py`

**Step 1: Write failing tests**

```python
import pandas as pd
import numpy as np
import pytest
from src.contact_analysis import extract_foldon_contacts, compute_all_foldon_cos

def test_extract_foldon_contacts_f01():
    """F01 segment 2 should have 22 L2 contacts (from segment_types)"""
    contacts_df = pd.read_csv('data/upstream/contacts_awsem.csv')
    seg_df = pd.read_csv('data/upstream/segment_types.csv')

    f01_contacts = contacts_df[contacts_df['family_id'] == 'F01']
    f01_segs = seg_df[seg_df['protein_id'] == 'F01']

    seg2 = f01_segs[f01_segs['seg_idx'] == 2].iloc[0]
    l2_pairs = extract_foldon_contacts(f01_contacts, seg_idx=2)
    assert len(l2_pairs) == seg2['n_L2']  # should be 22

def test_compute_all_foldon_cos_f01():
    """F01 should have CO for type A segments only (1, 2, 3), None for type B (0, 4)"""
    contacts_df = pd.read_csv('data/upstream/contacts_awsem.csv')
    seg_df = pd.read_csv('data/upstream/segment_types.csv')

    result = compute_all_foldon_cos('F01', contacts_df, seg_df)
    # result is list of dicts: [{seg_idx, seg_type, co, kf, n_l2, seg_start, seg_end, seg_size}, ...]
    assert len(result) == 5  # 5 segments

    # type B segments (0, 4) should have co=None, kf=None
    assert result[0]['co'] is None
    assert result[4]['co'] is None

    # type A segments (1, 2, 3) should have co > 0 and kf > 0
    for idx in [1, 2, 3]:
        assert result[idx]['co'] is not None
        assert result[idx]['co'] > 0
        assert result[idx]['kf'] > 0
```

**Step 2: Run to verify failure**

**Step 3: Write src/contact_analysis.py**

Functions:
- `extract_foldon_contacts(contacts_df, seg_idx)`:
  - Filter contacts_df where layer='L2' and seg_i == seg_j == seg_idx
  - Return list of (i, j) tuples

- `compute_all_foldon_cos(protein_id, contacts_df, seg_df)`:
  - For each segment of the protein:
    - Extract L2 contacts for that segment
    - Compute CO using contacts relative to seg_start (i.e., i_local = i - seg_start)
    - Compute kf from CO (Plaxco)
    - Return list of dicts with seg_idx, seg_type, co, kf, n_l2, seg_start, seg_end, seg_size

Important: contact residue indices i,j in contacts_awsem.csv are 0-indexed relative to the protein (not the foldon). For CO calculation, use the raw |i-j| values — the Plaxco formula uses absolute sequence separation, and L_i is foldon length.

**Step 4: Run tests, verify pass**

**Step 5: Commit**

```bash
git add src/contact_analysis.py tests/test_contact_analysis.py
git commit -m "feat: per-foldon contact order extraction from AWSEM contacts"
```

---

## Task 5: Rate Computation — Three Regimes

**Files:**
- Create: `cotrans-layer/src/rate_computation.py`
- Create: `cotrans-layer/tests/test_rate_computation.py`

**Step 1: Write failing tests**

```python
import numpy as np
import pytest
from src.rate_computation import uniform_rates, tai_rates, shuffled_synonymous_rates

def test_uniform_rates():
    """all codons get tau=0.06s"""
    codons = ['ATG', 'AAA', 'GCT', 'TGA']
    rates = uniform_rates(codons, tau_mean=0.06)
    assert len(rates) == 4
    assert np.allclose(rates, 0.06)

def test_tai_rates_mean():
    """mean of tAI rates should equal tau_mean"""
    codons = ['ATG', 'AAA', 'GCT', 'GAT', 'TTT'] * 20  # 100 codons
    tai_dict = {'ATG': 0.22, 'AAA': 1.0, 'GCT': 0.197, 'GAT': 0.295, 'TTT': 0.197}
    rates = tai_rates(codons, tai_dict, tau_mean=0.06, tai_floor=0.01)
    assert abs(np.mean(rates) - 0.06) < 0.001

def test_tai_rates_inverse_relationship():
    """higher tAI → lower tau (faster translation)"""
    codons = ['AAA', 'TTT']  # tAI: 1.0 vs 0.197
    tai_dict = {'AAA': 1.0, 'TTT': 0.197}
    rates = tai_rates(codons, tai_dict, tau_mean=0.06, tai_floor=0.01)
    assert rates[0] < rates[1]  # AAA faster

def test_shuffled_preserves_aa():
    """shuffled synonymous should change codon identity but preserve AA"""
    codons = ['GCT', 'GCC', 'GCA', 'GCG']  # all Ala
    tai_dict = {'GCT': 0.197, 'GCC': 0.333, 'GCA': 0.5, 'GCG': 0.16}
    rates = shuffled_synonymous_rates(codons, tai_dict, tau_mean=0.06,
                                       tai_floor=0.01, seed=42)
    # should still have 4 rates (one per position)
    assert len(rates) == 4
    # mean should still be ~0.06
    assert abs(np.mean(rates) - 0.06) < 0.01
```

**Step 2: Run to verify failure**

**Step 3: Write src/rate_computation.py**

Functions:
- `uniform_rates(codons, tau_mean=0.06)`:
  - Return np.full(len(codons), tau_mean)

- `tai_rates(codons, tai_dict, tau_mean=0.06, tai_floor=0.01)`:
  - For each codon, get tAI value (default to median tAI if codon not found)
  - Clip tAI to minimum of tai_floor (handles CGA=0.00007 → unreasonable tau)
  - Compute raw tau_i = 1.0 / tAI_i
  - Normalize: tau_i *= tau_mean / mean(tau_raw) so that mean tau = tau_mean
  - Return array

- `shuffled_synonymous_rates(codons, tai_dict, tau_mean=0.06, tai_floor=0.01, seed=None)`:
  - Group codons by amino acid (using CODON_TABLE from utils)
  - Within each AA group, randomly permute codon assignments (preserves AA identity)
  - Compute tAI rates on the shuffled codons
  - Return array

- `codon_table_lookup`: use CODON_TABLE from utils.py

**Step 4: Run tests, verify pass**

**Step 5: Commit**

```bash
git add src/rate_computation.py tests/test_rate_computation.py
git commit -m "feat: three rate regimes — uniform, tAI, shuffled-synonymous"
```

---

## Task 6: CDS Mapping — Port from Foldon Project

**Files:**
- Create: `cotrans-layer/src/cds_mapping.py`
- Create: `cotrans-layer/scripts/02_fetch_cds.py`

**Step 1: Write src/cds_mapping.py**

Port these functions from `/storage/kiran-stuff/foldon_project/codon_layer_analysis.py`:
- `fetch_sifts_mapping(pdb_id)` — adapted to use cotrans-layer cache paths
- `fetch_uniprot_cds_refs(accession)` — same
- `fetch_cds_from_ncbi(nucleotide_id, protein_id)` — same
- `fetch_cds_for_protein(pdb_id)` — same
- `align_cds_to_pdb(cds_protein, pdb_sequence)` — same
- `extract_codons(cds_sequence)` — new: split CDS into list of codon triplets

Key adaptations:
- Cache dirs: use both upstream symlinked caches AND local `data/cds/` cache
- Check upstream caches first (most SIFTS/CDS already cached in foldon_project)
- New CDS downloads go to `data/cds/{protein_id}_cds.json`
- Entrez.email = "research@cotrans-layer.org"

**Step 2: Write scripts/02_fetch_cds.py**

Reads data/protein_manifest.csv.
For each protein with cds_status='pending':
1. Call fetch_cds_for_protein(pdb_id)
2. If success: save CDS as FASTA to data/cds/{protein_id}.fasta
3. Update manifest with cds_status='success' or 'failed'
4. Log failures with reasons

Summary output: how many succeeded, how many failed, which failed.

**Step 3: Run it**

```bash
cd /storage/kiran-stuff/codon_project/cotrans-layer
python scripts/02_fetch_cds.py
```

Expected: ~26 from cache + some from new fetches. Target: ≥30/37 success.

**Step 4: Verify**

```bash
ls data/cds/*.fasta | wc -l
grep 'success' data/protein_manifest.csv | wc -l
```

**Step 5: Commit**

```bash
git add src/cds_mapping.py scripts/02_fetch_cds.py data/protein_manifest.csv
git commit -m "feat: CDS fetching pipeline, ported from foldon project"
```

---

## Task 7: Integration — Run Kinetic Model for All Proteins

**Files:**
- Create: `cotrans-layer/scripts/03_run_kinetic_model.py`

**Step 1: Write scripts/03_run_kinetic_model.py**

For each protein in manifest with cds_status='success':

1. Load segment_types for this protein
2. Compute per-foldon CO and kf (from Task 4)
3. Load CDS and extract codon list
4. For each of 36 parameter combinations (3 rates × 3 tunnels × 4 ku/kf):
   a. Compute codon_times array using the rate regime
   b. For each foldon: compute P_folded curve
   c. Record P_folded at self-emergence for each foldon
5. Save results:
   - `results/folding_curves/{protein_id}_cotrans.csv`: columns = chain_length, P_folded_seg_0, P_folded_seg_1, ... (for default params: uniform, tunnel=30, ku/kf=0.001)
   - `results/sensitivity/{protein_id}_sensitivity.csv`: columns = seg_idx, seg_type, co, kf, rate_regime, tunnel, ku_kf_ratio, p_folded_at_emergence, p_folded_final

For shuffled-synonymous: run 100 shuffles per protein, record mean and std of P_folded at emergence.

Use multiprocessing.Pool(workers=30) across proteins.

**Step 2: Run on 3 validation proteins first (checkpoint)**

```bash
python scripts/03_run_kinetic_model.py --proteins F38,F09,F01 --verbose
```

F38 = 1ubq (ubiquitin, 76 residues, simple), F09 = 8dfr (DHFR, 186 residues), F01 = 2abd (86 residues).

Check:
- P_folded ranges for type A foldons (expect 0.1-0.95)
- kf distribution (log-scale, expect 1-10000 s^-1 for small fragments)
- Curves look reasonable (monotonically increasing, plateau at P_eq)

**Step 3: If checkpoint passes, run all**

```bash
python scripts/03_run_kinetic_model.py --all --workers 30
```

**Step 4: Commit**

```bash
git add scripts/03_run_kinetic_model.py results/
git commit -m "feat: batch kinetic model — folding curves for all proteins"
```

---

## Task 8: L2/L3 Ordering Analysis

**Files:**
- Create: `cotrans-layer/scripts/04_analyze_layers.py`

**Step 1: Write scripts/04_analyze_layers.py**

Load all sensitivity CSVs from results/sensitivity/.

**Test 1: Type A Autonomous Folding Capacity**
- For each rate regime × tunnel × ku/kf combination:
  - Compute fraction of type A foldons with P_folded > 0.5 at self-emergence
  - Compute median P_folded at emergence
  - Report distribution (min, 25th, median, 75th, max)

**Test 2: L2-Before-L3 Temporal Ordering (core test)**
- For each protein, identify L3 interfaces: pairs of segments where contacts_awsem has L3 contacts between seg_i and seg_j (seg_i ≠ seg_j)
- For each L3 interface (seg_i, seg_j):
  - Self-folding time: each segment's P_folded at its own emergence
  - Interface time: chain_length = max(seg_end_i, seg_end_j) + tunnel
  - P_interface_possible = min(P_folded_i, P_folded_j) at interface_time
  - Compare: is P_folded(earlier_segment) > 0 when the interface becomes possible?
  - Core metric: fraction of L3 interfaces where BOTH participating foldons have P_folded > 0.5 at interface emergence (meaning L2 contacts already formed)
- Across proteins: Wilcoxon signed-rank test

**Test 3: Rate-Regime Invariance**
- For each foldon, compare P_folded(uniform) vs P_folded(tAI) vs P_folded(shuffled)
- Scatter plots: P_folded(uniform) vs P_folded(tAI) per foldon
- Pearson correlation (expect r > 0.95 if architecture-determined)
- Report: any foldons where ordering changes between regimes?

**Output files:**
- `results/summary/test1_autonomous_folding.csv`
- `results/summary/test2_l2_before_l3.csv`
- `results/summary/test3_regime_invariance.csv`
- `results/summary/analysis_report.md` — human-readable summary with all p-values and effect sizes

**Step 2: Run**

```bash
python scripts/04_analyze_layers.py
```

**Step 3: Review results/summary/analysis_report.md**

**Step 4: Commit**

```bash
git add scripts/04_analyze_layers.py results/summary/
git commit -m "feat: L2/L3 ordering analysis — all three tests"
```

---

## Task 9: Generate Figures

**Files:**
- Create: `cotrans-layer/scripts/05_generate_figures.py`

**Step 1: Write scripts/05_generate_figures.py**

Figure 1: **P_folded at emergence distribution**
- Three panels (one per rate regime)
- Within each: violin/strip plot of P_folded at emergence for all type A foldons
- Horizontal line at P_folded = 0.5
- Annotate with fraction > 0.5

Figure 2: **Co-translational folding curves — example proteins**
- 3-4 proteins selected for visual clarity (different sizes, different n_foldons)
- x = chain length (residues), y = P_folded
- Each foldon as a separate curve: type A in blue shades, type B as dashed red (flat at 0)
- Vertical dashed lines at foldon emergence points (seg_end + tunnel)
- Show default parameters (uniform, tunnel=30, ku/kf=0.001)

Figure 3: **Rate regime invariance scatter**
- x = P_folded(uniform), y = P_folded(tAI), one point per type A foldon
- Diagonal line for reference
- Pearson r annotated
- Color by protein

Figure 4: **Sensitivity heatmap**
- Rows: 36 parameter combinations
- Columns: metric (fraction_A_above_0.5, fraction_L3_ordered, pearson_uniform_tAI)
- Color: green (robust) to red (sensitive)

All figures saved to results/figures/ in both PNG (300 dpi) and PDF.

**Step 2: Run**

```bash
python scripts/05_generate_figures.py
```

**Step 3: Commit**

```bash
git add scripts/05_generate_figures.py results/figures/
git commit -m "feat: publication figures for cotrans-layer analysis"
```

---

## Task 10: Diagnostic Check and Final Report

**Files:**
- Modify: `cotrans-layer/results/summary/analysis_report.md`

**Step 1: Run diagnostic**

Check if the test is trivial:
- If >80% of type A foldons have P_folded > 0.95 at emergence under uniform rates:
  the test lacks discriminating power
- If >50% of type A foldons have P_folded < 0.1: model parameters are wrong
- Report the kf distribution (histogram) to show whether the Plaxco extrapolation
  to small foldons gives reasonable kinetics

**Step 2: Write final analysis_report.md**

Include:
- Dataset summary (N proteins, N foldons, N type A/B)
- Test 1 results: P_folded distribution, fraction > 0.5
- Test 2 results: L2-before-L3 ordering fraction, p-values, effect sizes
- Test 3 results: regime invariance (correlation, any differences?)
- Sensitivity: robust across 36 combinations?
- Diagnostic: is the test trivial or informative?
- Conclusion: one paragraph

**Step 3: Final commit**

```bash
git add -A
git commit -m "feat: complete cotrans-layer analysis pipeline"
```

---

## Execution Notes

### Dependencies between tasks

- Tasks 1-2: infrastructure, no dependencies
- Task 3: kinetic model core (pure math, no data needed)
- Task 4: contact analysis (needs upstream data symlinks from Task 1)
- Task 5: rate computation (needs utils from Task 1)
- Task 6: CDS mapping (needs manifest from Task 2, independent of Tasks 3-5)
- Task 7: integration (needs Tasks 3, 4, 5, 6 all complete)
- Task 8: analysis (needs Task 7)
- Task 9: figures (needs Task 8)
- Task 10: final (needs Tasks 8-9)

### Parallelizable tasks

Tasks 3, 4, 5, 6 are independent and can be developed in parallel after Tasks 1-2.

### Critical checkpoint

After Task 7, Step 2 (3-protein validation run): STOP and inspect results before proceeding. If P_folded values are unreasonable, diagnose before running all 37 proteins.

### Key file paths

| What | Path |
|------|------|
| Working dir | `/storage/kiran-stuff/codon_project/cotrans-layer/` |
| Upstream segments | `/storage/kiran-stuff/foldon_project/layer_results/segment_types.csv` |
| Upstream contacts | `/storage/kiran-stuff/foldon_project/layer_results/contacts_awsem.csv` |
| SIFTS cache | `/storage/kiran-stuff/foldon_project/uniprot_cache/` |
| CDS cache | `/storage/kiran-stuff/foldon_project/codon_results/cds_cache/` |
| tAI data | `/storage/kiran-stuff/codon_project/ecoli_tai_ws.tsv` |
| mu data | `/storage/kiran-stuff/codon_project/codon_error_rates.tsv` |
