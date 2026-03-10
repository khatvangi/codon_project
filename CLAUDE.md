# Codon Project

## Overview

This project investigates **codon-level biology** — how organisms select synonymous codons and what functional consequences that selection has. It builds on two prior projects:

1. **Proteostasis Law** (`/storage/kiran-stuff/proteostasis_law/`) — established the (mu, tAI) codon framework
2. **Foldon Project** (`/storage/kiran-stuff/foldon_project/`) — applied codon analysis to protein structural layers (result: negative)

The codon project is the NEXT phase — new questions beyond what was tested in foldon.

---

## GUIDELINES

Same as parent: `/storage/kiran-stuff/.claude/CLAUDE.md`

---

## Key Concepts

### The (mu, tAI) Framework
Every codon has two measurable properties:
- **mu** (mistranslation rate): probability of amino acid misincorporation per codon per translation event. From Landerer et al. (MBE). Range: 3.3e-05 (GAT) to 0.02 (CCC)
- **tAI** (tRNA Adaptation Index): decoding speed/efficiency based on tRNA gene copy numbers. Range: ~0.0001 (CGA) to 1.0 (AAA)

These define 4 **codon modes** per amino acid (using per-AA median thresholds):
- **Q1**: safe_sprinter — low mu, high tAI (accurate AND fast)
- **Q2**: safe_careful — low mu, low tAI (accurate but slow)
- **Q3**: risky_sprinter — high mu, high tAI (fast but error-prone)
- **Q4**: risky_careful — high mu, low tAI (worst of both)

### RSCU (Relative Synonymous Codon Usage)
RSCU = observed_count / expected_count for synonymous codons. RSCU=1 means no bias. RSCU>1 means preferred. WARNING: raw RSCU comparisons across protein regions are confounded by amino acid composition (hydrophobic AAs have different baseline RSCU).

### Deconfounding
- **Relative adaptiveness**: w = RSCU / RSCU_max (within-AA normalized, removes composition bias)
- **Stratified per-AA test**: compare metric for same amino acid at different positions, combine via sign test
- **mu_normalized**: mu(codon) / mean_mu(amino_acid) — removes AA identity from accuracy metric
- **Same-AA mode switching**: the zero-confound test — compare codons only for the same amino acid at different structural positions

---

## What Was Already Done (in Foldon Project)

### Test 1: Exon boundaries vs foldon boundaries — NEGATIVE
- 37 Galpern proteins, 156 exon boundaries vs AWSEM layer transitions
- Mean distance to nearest boundary: completely random (p uniformly distributed)
- Exon structure and protein folding modules are unrelated

### Test 2: RSCU by structural layer — PARTIALLY CONFOUNDED
- Raw: L2 RSCU > L3 RSCU (p=0.0003) — BUT ~half was AA composition artifact
- Deconfounded (relative adaptiveness w): L2=0.826 > L3=0.803 (p=0.025, survives)
- Per-AA sign test: 15/18 AAs show L2 > L3 (83%, p=0.008)
- Conclusion: weak but real signal — L2 positions use slightly more "preferred" codons

### Test 3: RSCU boundary pause — NEGATIVE
- No RSCU dip at foldon boundaries (p=0.45)

### Test 4: GC content at boundaries — MARGINAL
- p=0.026, not convincing

### Test 5: mu by layer (deconfounded) — NEGATIVE
- L2 vs L3 deconfounded mu: p=0.69 (coin flip)
- Per-protein: 11/24 (46%) show predicted direction

### Test 6: tAI by layer — BORDERLINE
- Deconfounded tAI: L2=1.117 > L3=1.072 (p=0.12)
- Per-protein: 17/24 (71%), p=0.064

### Test 7: Quadrant distribution by layer — NEGATIVE
- Chi-square p=0.12

### Test 8: Same-AA mode switching (THE KILLER TEST) — NEGATIVE
- mu: 8/18 AAs (44%) show L2 < L3, Fisher p=0.43
- tAI: 11/18 (61%), Fisher p=0.48
- **No layer-aware codon selection exists** for these 24 proteins

### Test 9: Error budget by layer — OPPOSITE DIRECTION
- j_L4 > j_L2 > j_L3 > j_surface (opposite to prediction)
- L4 inflated by Cys enrichment at disulfide bonds (Cys has highest mu)

### Bottom Line
The foldon structural layers do NOT drive codon selection. There is a weak RSCU/tAI signal that survives deconfounding, but the definitive same-AA mode switching test is flat. The codon axis is **closed for structural layers**.

---

## Data Files

### In this directory (to be populated)
New analyses go here.

### Source data (from prior projects)
| File | Location | Description |
|------|----------|-------------|
| `codon_error_rates.tsv` | `foldon_project/` or `proteostasis_law/errors/` | Per-codon mu (61 codons) |
| `ecoli_tai_ws.tsv` | `foldon_project/` or `proteostasis_law/errors/` | Per-codon tAI (61 codons) |
| `aa_mode_summary.tsv` | `proteostasis_law/errors/` | Per-AA mode counts |
| `codon_modes_ecoli.tsv` | `proteostasis_law/errors/` | E. coli codon mode assignments |
| `global_codon_usage.tsv` | `proteostasis_law/errors/` | Global codon usage frequencies |
| `residue_kappa_table.tsv` | `proteostasis_law/errors/` | Per-residue kappa values |
| `uniprot_pos_codon_ecoli.tsv` | `proteostasis_law/errors/` | E. coli proteome codon positions |
| `residue_codon_table.csv` | `foldon_project/codon_mode_results/` | 3976 residues with codon, layer, mu, tAI, quadrant |
| `contacts_awsem.csv` | `foldon_project/layer_results/` | 38,027 AWSEM contacts with L2/L3 classification |

### Key external databases
- CDS sequences: PDB → SIFTS → UniProt → EMBL → NCBI GenBank
- CDS cache: `foldon_project/codon_results/cds_cache/` (30 proteins already cached)

---

## Prior Project Scripts (reference, don't modify)

### Foldon project codon scripts
| Script | What it does |
|--------|-------------|
| `foldon_project/codon_layer_analysis.py` | RSCU analysis, CDS fetching pipeline, 3 signal tests |
| `foldon_project/codon_rscu_deconfound.py` | AA composition deconfounding (w, stratified tests) |
| `foldon_project/codon_mode_analysis.py` | Full (mu, tAI) framework, 7 analysis steps, mode switching |
| `foldon_project/codon_layer_test1.py` | Exon boundary vs foldon boundary null test |

### Proteostasis law scripts
| Script | What it does |
|--------|-------------|
| `proteostasis_law/figure4_codon_modes.py` | Original codon mode visualization |
| `proteostasis_law/figure4_corrected_null.py` | Corrected null model for mode distribution |
| `proteostasis_law/figure5_mode_richness.py` | Mode richness analysis |
| `proteostasis_law/figure6_capacity_bound.py` | Theoretical capacity bounds |
| `proteostasis_law/verify_diversity_claims.py` | Verification of diversity statistics |

---

## Protein Layer Architecture (from Foldon Project)

For reference — these are the structural layers that the codon analysis was tested against:

- **L2 (intra-foldon)**: contacts within a folding segment. E_direct < 0 (favorable). More conserved (p=0.005). These are the "autonomous core" contacts.
- **L3 (inter-foldon)**: contacts between folding segments. E_direct > 0 (unfavorable on average). Less conserved. These are the "cooperative interface" contacts.
- **L4 (functional)**: contacts at catalytic/binding sites. Not a separate energetic tier but a separate conservation tier (most constrained on both identity and class entropy).
- **Surface**: residues not in any contact layer.

Key finding: L2 and L3 form a bimodal energy distribution (ΔBIC=36,103). The boundary is a step function (transition width = 0.3 residues). 96% of proteins (105/109 tested) show L2 < L3 energy.

---

## Environment

- Machine: boron (10.147.17.21)
- Conda: use `openawsem` env if AWSEM tools needed, otherwise base
- Python: standard scientific stack (numpy, scipy, pandas, matplotlib, biopython)
- SLURM: use for jobs > 5 min
- 64 CPUs + 2 GPUs available

---

## Circular Permutation Analysis (ACTIVE)

Tests whether the DP-optimal modular decomposition is an emergent consequence of which residues are contiguous in the linear chain. Circular permutation changes contiguity without changing 3D structure.

### Implementation
- Script: `cp_analysis.py` in this directory
- **Permutation trick**: permute AWsemData arrays to CP coordinates, recompute pair_mask for CP backbone connectivity, run standard LayerScan pipeline unmodified
- AWSEM project for 1RIS: `foldon_project/cp_analysis/1ris_A/1ris/`
- AWsemData cache: `foldon_project/cp_analysis/1ris_results/1ris.npz`
- WT pickle: `foldon_project/cp_analysis/1ris_results/1ris/wt_result.pkl`
- CP results output: `cp_results/` in this directory

### S6 (PDB 1RIS) WT Results
- 97 residues, 5 foldons (F0-F4), 28 L2 / 447 L3 contacts
- F2 [44-71] is the strongest autonomous module: 17 L2 contacts, E_direct = -0.186

### CP54 Results (Lindberg et al. variant)
- Cuts through WT-F2, the strongest module
- 6 foldons (vs 5 in WT), 11 L2 / 466 L3 contacts
- WT-F2 destroyed: split into CP-F0 (WT[54-79]) and CP-F5 (WT[43-53])
- Novel wrapping module: CP-F2 = WT[0-8 + 92-96]
- 77/496 contacts switched layers (15.5%)
- L2→L3 contacts carry anomalously favorable E_direct (-0.119 vs typical L3 +0.029)

### S6 CP Variants
| Variant | Cut Site | Status | Notes |
|---------|----------|--------|-------|
| CP54 | 54 | DONE | strongest phenotype, cuts WT-F2 |
| CP13 | 13 | DONE | control, cuts in WT-F1 |
| CP68 | 68 | DONE | control, cuts near end of WT-F2 |

### S6 Insulation Analysis
- WT insulation: 0.995 (nearly perfect factorizability)
- CP54 insulation: **0.929** (14× more coupling — factorizability broken)
- CP13: 0.969, CP68: 0.978 (closer to WT)
- WT weak segments = the two Type A foldons (F1 ρ=0.105, F2 ρ=-0.129)
- CP54 destroys F2; F1 core → CP-F3 with amplified coupling (ρ=0.205)
- Novel wrapping module CP-F2 also weakly insulated (ρ=0.153)
- Key insight: gap widens (DP concentrates best contacts) BUT insulation drops (misplaced contacts create cross-layer correlations)

### S6 Figures
- Main text: `cp_results/cp_main_figure.png/pdf` — 3-column (arcs + mismatch bars + insulation)
- SI: `cp_results/cp_all_variants_si.png/pdf` — 4-row (WT, CP13, CP54, CP68)

### T4 Lysozyme CP Analysis (2LZM)
Second protein for CP analysis. 164 residues, two-subdomain architecture.
AWSEM project: `foldon_project/cp_analysis/2lzm_A/2lzm/`
AWsemData cache: `foldon_project/cp_analysis/2lzm_A/2lzm_awsem.npz`
WT result pickle: `foldon_project/cp_analysis/2lzm_A/wt_result.pkl`
Script: `cp_figure_t4l.py`, output: `cp_results_t4l/`

### T4L WT Results
- 164 residues, 9 foldons (4 Type A, 5 Type B), 113 L2 / 614 L3 contacts
- F0 [0-37]: helix A + N-sub N-terminus (38 res, 63 L2, E=+0.039)
- F2 [48-75]: N-subdomain core (28 res, 11 L2, E=-0.122)
- F6 [106-131]: C-subdomain core (26 res, 20 L2, E=-0.135, strongest)
- F8 [142-163]: C-terminal segment (22 res, 19 L2, E=+0.023)
- Insulation: 0.991

### T4L CP Results
| Variant | Exp. ΔΔG | Foldons | L2 | Switched | Insulation |
|---------|----------|---------|-----|----------|------------|
| WT | 0 | 9 | 113 | — | 0.991 |
| CP37 | 0.8 kcal/mol | 9 | 53 | 157 (20.9%) | 0.987 |
| CP13 | 3.0 kcal/mol | 9 | 84 | 100 (13.6%) | 0.959 |
| CP75 | 9.0 kcal/mol | 9 | 105 | 81 (10.8%) | 0.927 |

Key findings:
- Insulation perfectly predicts experimental stability loss: r=-0.980, p=0.020
- CP37 (mildest): most switching but highest insulation — clean reorganization at boundary
- CP75 (severest): least switching but lowest insulation — breaks module coupling
- CP13 wrapping module (F8' = WT[0-12 + 159-163]) recombines helix A with C-terminal tail, matching structural biology (helix A structurally belongs to C-domain)

### T4L Figures
- Main text: `cp_results_t4l/t4l_cp_main_figure.png/pdf` — 3-column (arcs + insulation bars + insulation vs ΔΔG scatter)
- SI: `cp_results_t4l/t4l_cp_all_variants_si.png/pdf` — 4-row (WT, CP37, CP13, CP75)

### T4L Literature
- Zhang, Baase & Matthews 1993 (Biochemistry) — CP37, 0.8 kcal/mol loss
- Llinas & Marqusee 1998 (Protein Science) — CP13 (3 kcal/mol), CP75 (9 kcal/mol)
- Cellitti et al. 2007 (Protein Science) — CP13* crystal structure (PDB 2O4W)
- Shank et al. 2010 (Nature) — CP13 abolishes mechanical cooperativity

### Key Technical Notes
- `pair_mask` must be recomputed for CP (encodes polymer connectivity, not 3D distance)
- Burial weights recomputed via `precompute_contact_weights` on permuted coords
- "absent" contacts = pairs whose pair_mask status changed at cut/junction sites
- STRIDE not needed for contact energy analysis (only for MD ssweight)
- gamma.dat and burial_gamma.dat are protein-independent AWSEM parameters

---

## Open Questions (potential directions)

These are starting points — the user will define the actual research direction:

1. **Cross-species codon selection**: do the (mu, tAI) patterns found in E. coli generalize to other organisms? The proteostasis_law project has `cross_species/` data.
2. **Codon mode richness**: what determines how many modes an amino acid uses across the proteome? (started in proteostasis_law)
3. **Capacity bounds**: theoretical limits on error correction via codon selection (started in proteostasis_law)
4. **Codon co-occurrence / mRNA structure**: do adjacent codons show correlated selection (independent of protein structure)?
5. **Codon selection at protein-protein interfaces vs cores**: different structural question from L2/L3 (which are folding layers, not interface layers)
6. **Translation speed and cotranslational folding**: does tAI modulate folding kinetics at domain boundaries?
