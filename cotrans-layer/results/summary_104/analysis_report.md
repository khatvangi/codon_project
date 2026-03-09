# emergence gap analysis — extended to 104 proteins

## overview

extension of the 37-protein Galpern analysis to the full 104-protein
dataset from table S1. protein_L (2ptl) excluded: all type B, no contacts
in energy_decomposition.csv.

## dataset summary

- **proteins**: 104 (37 Galpern + 67 validation)
- **L3 interfaces**: 3071
- **type A foldons with L3 interfaces**: 400
- **type B foldons**: 549

## comparison: 37 Galpern vs 104 full set

| metric | 37 Galpern | 104 full | delta |
|--------|-----------|----------|-------|
| L3 interfaces | 1244 | 3071 | +1827 |
| median gap (residues) | 54 | 48 | -6 |
| min gap (residues) | 10 | 10 | +0 |
| mean gap (residues) | 78.0 | 65.3 | -12.6 |
| zero-gap interfaces | 0/1244 (0.0%) | 0/3071 (0.0%) | |
| type A within empirical kf | 165/165 (100.0%) | 400/400 (100.0%) | |
| temporal order OK | 37/37 (100.0%) | 104/104 (100.0%) | |

## emergence gap distribution (tunnel = 30)

| statistic | residues | seconds |
|-----------|----------|---------|
| min       | 10 | 0.600 |
| 25th      | 24 | 1.440 |
| median    | 48 | 2.880 |
| 75th      | 87 | 5.220 |
| max       | 360 | 21.600 |
| mean      | 65.3 | 3.921 |

- **zero-gap interfaces**: 0/3071 (0.0%)

## minimum kf required (tunnel = 30)

| statistic | log10(kf_min) |
|-----------|---------------|
| min       | -0.72 |
| 25th      | -0.19 |
| median    | -0.08 |
| 75th      | -0.02 |
| max       | 0.06 |

- **within empirical range**: 400/400 (100.0%)

## temporal classification (tunnel = 30)

- **proteins with valid temporal ordering**: 104/104 (100.0%)
  all proteins pass.

## sensitivity to tunnel length

| tunnel | median gap (res) | median gap (s) | foldons within kf | temporal OK |
|--------|-----------------|----------------|-------------------|-------------|
| 25 | 48 | 2.880 | 400/400 (100.0%) | 104/104 (100.0%) |
| 30 | 48 | 2.880 | 400/400 (100.0%) | 104/104 (100.0%) |
| 40 | 48 | 2.880 | 400/400 (100.0%) | 104/104 (100.0%) |

## key finding

the 37-protein result generalizes to the full 104-protein dataset.
the emergence gap is a robust, parameter-free geometric property of
protein modular architecture. the layer architecture is a construction
plan compatible with the ribosome's vectorial N→C synthesis.
