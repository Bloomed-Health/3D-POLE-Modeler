# 3D POLE Modeler

[![CI](https://github.com/bloomed-health/3D-POLE-Render/actions/workflows/ci.yml/badge.svg)](https://github.com/bloomed-health/3D-POLE-Render/actions/workflows/ci.yml)
[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](LICENSE)
[![Node](https://img.shields.io/badge/node-%3E%3D20-brightgreen.svg)](https://nodejs.org/)
[![Contributing](https://img.shields.io/badge/PRs-welcome-brightgreen.svg)](CONTRIBUTING.md)

Interactive 3D structural viewer for the human DNA Polymerase epsilon (Pol e) catalytic subunit (p261, UniProt Q07864).

## Overview

This tool renders a **procedural cartoon-ribbon schematic** of the POLE holoenzyme directly in the browser using Three.js (WebGL). It is designed for exploring the domain architecture, catalytic mechanism, cancer-associated mutations, and clinical consequences of rare germline variants of Pol e.

> **Important:** This is a schematic teaching model, not an atomic-coordinate viewer. The ribbon geometry, domain centers, and secondary-structure elements are procedurally generated to illustrate domain topology and functional relationships. Color modes (B-factor, conservation, charge, pathogenicity) display approximate domain-level or synthetic per-residue values, not crystallographic or computed data. For analysis of real atomic coordinates, load PDB [9F6D](https://www.rcsb.org/structure/9F6D) (Roske & Yeeles 2024, human Pol e-PCNA-DNA) into [Mol\*](https://molstar.org/) or PyMOL.

## Features

- **Cartoon ribbon representation** - helices, beta-strands (with arrowheads), and loops rendered per secondary structure element across all eight functional domains (NTD, EXO, Palm, P, Fingers, Thumb, Inactive-Pol, CTD)
- **Two-metal-ion catalysis** - ball-and-stick active site with octahedral Mg2+ coordination, conserved aspartate residues (D640, D642, D860), and distance annotations
- **Incoming dNTP + Watson-Crick base pairing** - template nucleotide with hydrogen bond visualization
- **3'-to-5' exonuclease site** - two-metal proofreading active site with catalytic residues D275, E277, D368, D370
- **Zinc-finger motifs** - CysA/CysB Zn2+ tetrahedral coordination in the CTD
- **[4Fe-4S] cluster** - iron-sulfur cluster in the P domain (CysX motif), a feature unique to Pol e among B-family polymerases
- **DNA double helix** - template and primer backbone tubes with base-pair rungs
- **Mutation hotspot markers** - P286R, V411L, S297F, L424V, D287E, P436R, M444K, S459F, F367S with interactive tooltips
- **Accessory subunits** - ghost representations of POLE2 (p59), POLE3 (p17), POLE4 (p12)
- **PCNA sliding clamp** - toroidal processivity factor with homotrimer symmetry markers and tripartite Pol e contact interface (PIP-box Q1180, Thumb insertion res. 1102-1122, P-domain contact)
- **Beta-sheet hydrogen bonds** - inter-strand H-bond visualization (schematic, not from DSSP)
- **Exo-Pol shuttling path** - animated primer terminus transfer between active sites (~40 A, per Roske & Yeeles 2024)

### Color Modes

| Mode | Description |
|------|-------------|
| Domain | Default domain-based coloring |
| B-factor | Approximate flexibility gradient by SS type (not crystallographic B-factors) |
| Conservation | Domain-level evolutionary conservation estimate (not per-residue MSA/ConSurf) |
| Variant | COSMIC/ClinVar variant density (domain-level approximation) |
| Charge | Amino acid formal charge (domain-level average, not APBS electrostatics) |
| Pathogenicity | Clinical significance spectrum (blue-green = benign, white = VUS, red = pathogenic) |

### Color Schemes

- **Default** - muted, publication-quality palette
- **PyMOL** - classic PyMOL-style bright colors
- **Accessible** - colorblind-friendly palette

## Mol* Coordinate Viewer

A separate page (`molstar.html`) provides a **real coordinate viewer** using [Mol*](https://molstar.org/), the same viewer used by RCSB PDB and PDBe. Features:

- Loads PDB 9F6D (Roske & Yeeles 2024) with domain-colored cartoon representation
- **Finger conformational triad**: switch between open (9F6D), ajar (9F6E), and closed (9F6F) states
- **PDB 9B8S** (He et al. 2024) — independent Pol e-PCNA cryo-EM structure
- Real DSSP secondary structure, crystallographic B-factors, metal coordination geometry
- Multiple representations: cartoon, ball-and-stick, surface, putty (B-factor), spacefill
- Domain coloring, chain coloring, element coloring, secondary structure coloring
- Toggle Mol* built-in controls for advanced operations

## Variant Catalogue

A variant analysis page (`variant.html`) provides per-mutation analysis for all characterized POLE pathogenic variants:

- **ACMG classification** with Bayesian posterior probabilities (Mur et al. 2023 framework)
- **Auditable evidence codes** — each PS/PM/PP code has a literature citation
- **Mutational signature attribution** — per-variant SBS10a/b/28/14 weights with bar charts
- **DDG stability prediction** — heuristic destabilization estimate (not FoldX/Rosetta)
- **Mol* structural context** — variant residue highlighted in PDB 9F6D coordinates
- URL-addressable variants: `variant.html#P286R`, `variant.html#V411L`, etc.

## Mutation Viewer

A separate page (`mutation.html`) visualizes the ultra-rare **c.138del** germline variant (p.Leu46Phefs*8), showing:

- Side-by-side wild-type vs. truncated mutant comparison
- Domain retention bar chart (54 / 2,286 residues = 2.4% translated)
- Reading frame analysis with frameshift coloring

**Clinical note:** truncating POLE variants such as c.138del are typically classified as likely benign for polymerase-proofreading-associated polyposis (PPAP). PPAP requires an active but proofreading-deficient polymerase to generate a mutator phenotype; a null allele (as produced by early frameshift) does not confer a mutator effect, and the wild-type allele provides sufficient function. See Mur et al. 2023 for the gene-specific ACMG/AMP classification framework.

## Julia Backend Pipeline

A Julia data pipeline (`pipeline/`) computes variant scoring, mutational signature deconvolution, Bayesian variant classification, and structural math from PDB coordinates, outputting JSON consumed by the viewer at runtime. See `pipeline/README.md` for details. When pipeline JSON is unavailable, the viewer falls back to built-in default values.

## Running Locally

Serve the project directory with any static HTTP server:

```bash
npx serve .
```

Then open `http://localhost:3000` in a modern browser.

## Controls

| Input | Action |
|-------|--------|
| Drag | Orbit camera |
| Shift+Drag | Pan |
| Scroll | Zoom |
| Hover | Inspect elements (tooltip) |
| `1`-`8` | Focus domain |
| `R` | Toggle auto-rotate |
| `S` | Screenshot (PNG) |
| `M` | Measure distance |
| `H` | Toggle helix H-bonds |
| `B` | Toggle beta-sheet H-bonds |
| `Z` | Toggle Zn-finger motifs |
| `P` | Toggle PCNA clamp |
| `D` | Toggle DNA |
| `L` | Toggle labels |
| `E` | Toggle Exo-Pol path |
| `?` | Keyboard shortcuts |

## Limitations

This viewer is a **schematic teaching tool**, not a molecular graphics program:

- **Coordinates are procedural.** Ribbon paths are generated algorithmically to approximate domain topology, not parsed from PDB ATOM records. Domain centers are hand-placed to reflect the spatial arrangement in the Roske & Yeeles 2024 structures but are not superimposed on deposited coordinates.
- **B-factors are not crystallographic.** The B-factor color mode assigns a single pseudo-value per secondary-structure type (helix = 0.2, strand = 0.5, loop = 0.9). Real B-factors require a refined coordinate set.
- **Conservation is domain-level.** Values are literature-derived per-domain estimates, not per-residue Shannon entropy from a multiple sequence alignment. For real per-residue conservation, use [ConSurf](https://consurf.tau.ac.il/) with UniProt Q07864.
- **Electrostatic surface is not computed.** The charge mode shows formal amino acid charge averaged per domain, not a Poisson-Boltzmann (APBS/PDB2PQR) electrostatic potential mapped onto a solvent-excluded surface.
- **The molecular envelope is a geometric approximation**, not a Connolly/MSMS solvent-excluded surface.
- **H-bonds are schematic.** They are placed along the backbone to illustrate i->i+4 and inter-strand patterns, but donor-acceptor distances and angles are not computed from atomic coordinates via DSSP.

For coordinate-derived molecular visualization, load PDB **9F6D** (human Pol e-PCNA-DNA, open Finger conformation) into [Mol\*](https://molstar.org/) or [PyMOL](https://pymol.org/).

## References

- UniProt [Q07864](https://www.uniprot.org/uniprot/Q07864) (POLE_HUMAN)
- PDB: [9F6D](https://www.rcsb.org/structure/9F6D), [9F6E](https://www.rcsb.org/structure/9F6E), [9F6F](https://www.rcsb.org/structure/9F6F) (Roske & Yeeles 2024); [9B8S](https://www.rcsb.org/structure/9B8S) (He et al. 2024); [4M8O](https://www.rcsb.org/structure/4M8O) (Hogg et al. 2014)
- Roske & Yeeles, *Nat Struct Mol Biol* 31:1921-1932 (2024) — first human Pol e catalytic domain structures: open (9F6D), ajar (9F6E), and closed (9F6F) Finger conformations; four mismatch proofreading intermediates (post-insertion, arrest, frayed substrate, mismatch excision)
- He, Wang, Yao, O'Donnell & Li, *Nat Commun* 15:7847 (2024) — independent human Pol e-PCNA cryo-EM structure (9B8S)
- Mur, Viana-Errasti et al., *Genome Med* 15:85 (2023) — gene-specific ACMG/AMP framework for POLE/POLD1 germline variant classification with likelihood-ratio-based Bayesian update
- Yuan, Georgescu, Schauer, O'Donnell & Li, *Nat Commun* 11:3156 (2020) — yeast Pol e holoenzyme structure revealing the non-catalytic Pol2 module
- Robinson, Coorens et al., *Nat Genet* 53:1434 (2021) — normal-tissue somatic mutation rates in POLE/POLD1 germline carriers (Sanger)
- Hogg, Osterman, Goldsmith et al., *Nat Struct Mol Biol* 21:49-55 (2014) — S. cerevisiae Pol2 catalytic domain (4M8O)

## License

[Apache 2.0](LICENSE) — see [SECURITY.md](SECURITY.md) for vulnerability reporting.
