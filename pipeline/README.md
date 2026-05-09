# POLE Pipeline

Julia data pipeline for computing variant scoring, mutational signature deconvolution,
Bayesian variant classification, and structural math from PDB coordinates.

## Requirements

- Julia 1.10+
- Internet access (for PDB download)

## Setup

```bash
cd pipeline
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

## Running

```bash
# Full pipeline (downloads PDB, computes all outputs)
julia --project=. run.jl

# Without PDB download (uses synthetic/built-in data only)
julia --project=. run.jl --no-pdb
```

## Output

JSON files are written to `../data/`:

| File | Contents |
|------|----------|
| `pole_structure.json` | Domain geometry, SS elements, metal sites, backbone coords |
| `pole_scores.json` | Per-residue B-factor, conservation, variant density, pathogenicity, charge |
| `pole_mutations.json` | Mutation annotations, SBS signature attribution, ACMG classification with auditable citations |
| `pole_meta.json` | Pipeline version, PDB resolutions, data provenance, checksums |

## Data Sources

See `data/DATA_SOURCES.md` for instructions on downloading real variant databases
(ClinVar, COSMIC, AlphaMissense, ConSurf MSA, COSMIC SBS v3.4 profiles).

When data files are absent, the pipeline uses built-in fallback values.

## Modules

| Module | Purpose |
|--------|---------|
| `config.jl` | Domain ranges, PDB IDs, mutation definitions |
| `pdb.jl` | mmCIF parsing, coordinate extraction |
| `structure.jl` | DSSP-style SS assignment, B-factors, domain geometry |
| `conservation.jl` | MSA parsing, per-residue Shannon entropy |
| `variants.jl` | ClinVar/COSMIC/AlphaMissense processing |
| `signatures.jl` | NMF/NNLS signature deconvolution (SBS10a/b/28/14) |
| `classification.jl` | Bayesian ACMG classification (Mur et al. 2023) |
| `electrostatics.jl` | Per-residue formal charge |
| `ddg.jl` | Heuristic DDG stability estimation |
| `output.jl` | JSON serialization |

## Tests

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Limitations

- **DDG estimates** are heuristic (not FoldX/Rosetta). Use for qualitative assessment only.
- **Signature profiles** are approximations. Download real COSMIC SBS v3.4 for production.
- **Conservation** without an MSA file produces synthetic domain-level estimates.
- **ACMG evidence** is pre-curated with literature citations but should be reviewed by a clinical geneticist.
