# Co-Translational Folding × Layer Architecture: Design Document

Date: 2026-03-03

## Goal

Test whether the layer architecture (L2 intra-foldon / L3 inter-foldon) predicts
co-translational folding order: do L2 contacts form as each foldon emerges from
the ribosome, while L3 contacts form only after multiple foldons have emerged?

## Key Insight

Run the kinetic model under three rate regimes (uniform, tAI, shuffled-synonymous).
If L2-before-L3 ordering is identical across all three, the ordering is
architecture-determined, not codon-determined — connecting to the Table 3 negative
codon results from the main paper.

## Scope

- 37 Galpern proteins with full contact + segment data
- 354 segments total (165 type A with L2 contacts, 189 type B without)
- contacts_awsem.csv provides individual (i,j) pairs for per-foldon contact order
- Non-Galpern proteins (22) lack individual contact data; deferred to future extension

## Data Sources

All from /storage/kiran-stuff/foldon_project/:
- segment_types.csv: foldon boundaries, seg_start/seg_end, n_L2, seg_type (A/B)
- contacts_awsem.csv: 38,027 contacts with layer, seg_i, seg_j, family_id
- codon_error_rates.tsv + ecoli_tai_ws.tsv: per-codon mu and tAI (in codon_project/)
- cds_cache/: ~30 proteins already cached

## Architecture

```
protein_manifest.csv (37 proteins)
    → CDS fetching (PDB → SIFTS → UniProt → EMBL → NCBI)
    → Rate computation (3 regimes: uniform, tAI, shuffled-synonymous)
    → O'Brien kinetic model (per-foldon P_folded curves)
    → L2/L3 ordering tests + statistics
    → Figures
```

## Kinetic Model

Per-foldon contact order:
    CO(foldon_i) = (1 / (L_i * N_i)) * sum(|r1 - r2|)
    where sum is over L2 contacts within segment i only
    L_i = foldon length (residues), N_i = number of L2 contacts

Folding rate (Plaxco 1998):
    ln(kf) = -17.0 * CO + 8.0

Unfolding rate:
    ku = kf * ratio, where ratio in {0.0001, 0.001, 0.01, 0.1}

Translation time per codon:
    Uniform:  tau = 0.06s for all codons
    tAI:      tau_i = tau_mean / tAI(codon_i), normalized so mean = 0.06s
    Shuffled: tau from tAI after randomly swapping synonymous codons, 100x

Tunnel length: default 30 residues, sensitivity at 25 and 40

Folding probability at chain length L:
    emergence = seg_end + tunnel_length
    for L > emergence:
        t_avail = sum(tau_i) for i from emergence to L
        P_folded(L) = kf/(kf+ku) * (1 - exp(-(kf+ku) * t_avail))

Time-zero approximation: folding only possible once complete domain clears tunnel.
This is conservative (underestimates P_folded for large foldons whose N-terminal
portions have been outside for many codons).

## Type A vs Type B Foldons

Type A (n_L2 > 0): has autonomous intra-module contacts. CO and kf calculable.
Type B (n_L2 = 0): no autonomous core. No independent folding capacity.
    Predicted to acquire structure only through post-emergence interface formation.

Type B P_folded = 0 is definitional, not a test result. The genuine tests are below.

## Tests

### Test 1: Type A Autonomous Folding Capacity
Do type A foldons achieve substantial P_folded at self-emergence?
Not guaranteed — depends on CO, kf, time available, ku/kf.
Report: distribution of P_folded at emergence across 165 type A foldons.
Success: majority achieve P_folded > 0.5.

### Test 2: L2-Before-L3 Temporal Ordering (THE CORE TEST)
For each L3 interface between foldons i and j:
    interface_emergence = max(seg_end_i, seg_end_j) + tunnel
    P_interface_possible = min(P_folded_i, P_folded_j) at interface_emergence
Compare self-folding (at own emergence) vs interface availability (at later emergence).
Prediction: self-folding precedes interface formation.
Could fail if foldons are interleaved or fold slowly relative to neighbor emergence.

### Test 3: Rate-Regime Invariance
L2-before-L3 ordering identical across uniform, tAI, shuffled-synonymous?
If yes: ordering is architecture-determined, codons are irrelevant to it.

### Prediction (descriptive, not a test):
Type B foldons remain unstructured at self-emergence and acquire structure
only after flanking type A foldons fold and provide L3 interfaces.

## Sensitivity Analysis

36 parameter combinations per protein:
    3 rate regimes x 3 tunnel lengths x 4 ku/kf ratios

## Diagnostic Check

If >80% of type A foldons have P_folded > 0.95 at emergence under uniform rates,
the test lacks discriminating power. Report but note triviality.

## Statistics

- Wilcoxon signed-rank across proteins for P_folded distributions
- Cohen's d effect sizes
- Bootstrap 10,000 resamples for confidence intervals
- Multiple testing correction (9 comparisons: 3 regimes x 3 tunnels)

## Key Figures

1. P_folded at emergence: distribution across type A foldons, one panel per regime
2. Co-translational folding curves: 3-4 example proteins, type A blue, type B red
3. Regime comparison: scatter P_folded(uniform) vs P_folded(tAI), should be diagonal
4. Sensitivity heatmap: robustness of L2-before-L3 across 36 parameter combinations

## Directory Structure

cotrans-layer/
    data/
        protein_manifest.csv
        cds/
        upstream/ (symlinks to foldon_project files)
    src/
        cds_mapping.py
        rate_computation.py
        kinetic_model.py
        contact_analysis.py
        utils.py
    scripts/
        01_build_manifest.py
        02_fetch_cds.py
        03_run_kinetic_model.py
        04_analyze_layers.py
        05_generate_figures.py
    results/
        folding_curves/
        sensitivity/
        figures/
        summary/
    tests/
        test_kinetic_model.py

## Validation Checkpoint

Before full run, validate on CI2, DHFR_ecoli, barnase:
- Check P_folded ranges are reasonable (0.1-0.95 for most type A)
- Check kf distribution makes sense for foldon-sized fragments
- If trivial (all >0.99) or broken (all <0.01), diagnose before scaling
