# 3D POLE Modeler

[![CI](https://github.com/Bloomed-Health/3D-POLE-Modeler/actions/workflows/ci.yml/badge.svg)](https://github.com/Bloomed-Health/3D-POLE-Modeler/actions/workflows/ci.yml)
[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](LICENSE)
[![Node](https://img.shields.io/badge/node-%3E%3D20-brightgreen.svg)](https://nodejs.org/)
[![Contributing](https://img.shields.io/badge/PRs-welcome-brightgreen.svg)](CONTRIBUTING.md)

Interactive 3D structural viewer for the human DNA Polymerase epsilon (Pol e) catalytic subunit (p261, UniProt Q07864).

<img width="1438" height="821" alt="homescreen-demo" src="https://github.com/user-attachments/assets/b372734d-d3e4-48a8-bfdc-32bf4352f554" />

## Overview

This tool provides interactive 3D structural viewers for exploring the domain architecture, catalytic mechanism, cancer-associated mutations, and clinical consequences of rare germline variants of Pol e. The home page renders **real atomic coordinates** from PDB 9F6D using [Mol\*](https://molstar.org/) in ball-and-stick representation. A separate schematic page (`schematic.html`) renders a **procedural cartoon-ribbon model** using Three.js.

> **Note:** The schematic viewer is a teaching model, not an atomic-coordinate viewer. Its ribbon geometry, domain centers, and secondary-structure elements are procedurally generated to illustrate domain topology. For analysis of real atomic coordinates, use the home page Mol\* viewer or load PDB [9F6D](https://www.rcsb.org/structure/9F6D) into PyMOL.

## Schematic Viewer Features

The schematic viewer (`schematic.html`) includes:

- **Cartoon ribbon representation** - helices, beta-strands (with arrowheads), and loops rendered per secondary structure element across all eight functional domains (NTD, EXO, Palm, P, Fingers, Thumb, Inactive-Pol, CTD)
- **Two-metal-ion catalysis** - ball-and-stick active site with octahedral Mg2+ coordination, conserved aspartate residues (D626 [Motif A], D860/D862 [Motif C]; human POLE1, UniProt Q07864, verified from He et al. 2024, PDB 9B8S), and distance annotations
- **Incoming dNTP + Watson-Crick base pairing** - template nucleotide with hydrogen bond visualization
- **3'-to-5' exonuclease site** - two-metal proofreading active site with the DEDDy catalytic set: D275, E277 (ExoI), D368 (ExoII), D462 (ExoIII, verified from PDB [9F6L](https://www.rcsb.org/structure/9F6L), Roske & Yeeles 2024), Y458 (conserved tyrosine)
- **Zinc-finger motifs** - CysA/CysB in the CTD. CysA binds Zn2+ (required for PCNA binding). CysB is rendered as Zn2+ per the human p59-p261 crystal structure, though Netz et al. (2012) originally assigned a [4Fe-4S] cluster to CysB in B-family polymerases — the viewer depicts the Zn model; see the [4Fe-4S] entry for the CysX/P-domain cluster
- **[4Fe-4S] cluster** - iron-sulfur cluster in the P domain via the CysX motif (coordinating cysteines C651, C654, C663, C747 in human POLE1), a feature unique to Pol e among B-family polymerases (ter Beek et al. 2019; Jain et al. 2014)
- **DNA double helix** - template and primer backbone tubes with base-pair rungs
- **Mutation hotspot markers** - somatic EDMs (P286R, V411L, S297F, A456P, P436R, S459F per Leon-Castillo 2020) and germline PPAP variants (L424V, S297F, F367S, M444K, S459F per Palles 2013 / Mur 2023) with interactive tooltips. Note: L424V (Leu->Val) is the germline PPAP founder; somatic L424I (Leu->Ile) is a distinct variant. D287E is included as a somatic recurrence but is not in the Leon-Castillo pathogenic set. 10 variants total in the variant catalogue
- **Accessory subunits** - ghost representations of POLE2 (p59), POLE3 (p17), POLE4 (p12)
- **PCNA sliding clamp** - toroidal processivity factor with homotrimer symmetry markers and tripartite Pol e contact interface (PIP-box Q1180, Thumb insertion res. 1102-1122, P-domain contact)
- **Beta-sheet hydrogen bonds** - inter-strand H-bond visualization (schematic, not from DSSP)
- **Exo-Pol shuttling path** - animated primer terminus transfer between active sites (~40 A, per Roske & Yeeles 2024). **Note:** the current animation shows a straight-line path; Roske & Yeeles demonstrate a multi-step DNA handover (Thumb -> P domain [frayed substrate] -> back to Thumb [mismatch excision], with 1-bp backtracking), with alpha-helix 3 of the Exo domain ("wedge helix") physically blocking the direct path. The nascent strand is held throughout by an Anchor loop (~residues 970-980, K974/R975/R976). Replacing the linear tween with a four-state morph across the proofreading intermediates (9F6I->9F6J->9F6K->9F6L) is a planned upgrade

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

## Mol* Coordinate Viewer (Home Page)

The home page (`index.html`) is a **real coordinate viewer** using [Mol*](https://molstar.org/), the same viewer used by RCSB PDB and PDBe. Features:

- Loads PDB 9F6D (Roske & Yeeles 2024) with domain-colored cartoon representation
- **Finger conformational triad**: switch between open (9F6D), ajar (9F6E), and closed (9F6F) states. Four additional mismatch proofreading intermediates (9F6I-9F6L) are deposited but not yet wired into the selector. **Construct note:** 9F6I and 9F6J use an exonuclease-dead construct (D275A/E277A); active-site labels at the exo site correspond to alanines, not native aspartate/glutamate catalytic residues
- **PDB 9B8S** (He et al. 2024) — independent Pol e-PCNA cryo-EM structure
- Real DSSP secondary structure, crystallographic B-factors, metal coordination geometry
- Multiple representations: cartoon, ball-and-stick, surface, putty (B-factor), spacefill
- Domain coloring, chain coloring, element coloring, secondary structure coloring
- Toggle Mol* built-in controls for advanced operations

## Variant Catalogue

A variant analysis page (`variant.html`) provides per-mutation analysis for all characterized POLE pathogenic variants:

- **ACMG classification** with Bayesian posterior probabilities (Tavtigian et al. 2018 likelihood-ratio framework, applied per Mur et al. 2023)
- **Auditable evidence codes** — each PS/PM/PP code has a literature citation
- **Mutational signature attribution** — per-variant SBS10a/b/28/14 weights with bar charts
- **DDG stability prediction** — heuristic destabilization estimate (not FoldX/Rosetta)
- **Mol* structural context** — variant residue highlighted in PDB 9F6D coordinates
- URL-addressable variants: `variant.html#P286R`, `variant.html#V411L`, etc.

## Mutation Viewer

A separate page (`mutation.html`) visualizes the ultra-rare **c.138del** germline variant (p.Leu46Phefs*8), showing:

- Side-by-side wild-type vs. truncated mutant comparison
- Domain retention bar chart (stop codon at position 53; 52 / 2,286 residues = 2.3% translated)
- Reading frame analysis with frameshift coloring
- Hypothesis panel (H1–H5) for the observed genotype-phenotype anomaly

**Clinical note:** c.138del is classified as **Pathogenic** by the carrier's clinical lab. The carrier presents with a PPAP-consistent phenotype (polyposis, multiple primaries). While truncating POLE variants are not canonical PPAP drivers (PPAP requires an active but proofreading-deficient polymerase), the observed phenotype creates a genotype-phenotype anomaly under active investigation. Five candidate mechanisms are evaluated in the hypothesis panel. **Note:** the Mur et al. 2023 gene-specific ACMG/AMP framework covers germline variants *in the exonuclease domain* of POLE/POLD1; c.138del (NTD, LOF) falls outside its scope by both criteria. PVS1 is explicitly excluded from the framework. Biallelic POLE LOF variants cause recessive syndromes (FILS, IMAGe-I) distinct from PPAP.

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
- **Polymerase active-site residues.** D626 (Motif A) and D860/D862 (Motif C) are the human POLE1 (UniProt Q07864) catalytic aspartates, verified from He et al. 2024 (PDB 9B8S, Nat Commun 15:7847). The earlier D640/D642 annotation was an erroneous carry-over of *S. cerevisiae* Pol2 numbering — human residues at positions 640/642 are Leu/Pro, not Asp.

For coordinate-derived molecular visualization, load PDB **9F6D** (human Pol e-PCNA-DNA, open Finger conformation) into [Mol\*](https://molstar.org/) or [PyMOL](https://pymol.org/).

## References

- UniProt [Q07864](https://www.uniprot.org/uniprot/Q07864) (POLE_HUMAN)
- PDB: [9F6D](https://www.rcsb.org/structure/9F6D), [9F6E](https://www.rcsb.org/structure/9F6E), [9F6F](https://www.rcsb.org/structure/9F6F) (Roske & Yeeles 2024, matched DNA: open/ajar/closed Finger); [9F6I](https://www.rcsb.org/structure/9F6I), [9F6J](https://www.rcsb.org/structure/9F6J), [9F6K](https://www.rcsb.org/structure/9F6K), [9F6L](https://www.rcsb.org/structure/9F6L) (Roske & Yeeles 2024, mismatch proofreading intermediates); [9B8S](https://www.rcsb.org/structure/9B8S) (He et al. 2024); [4M8O](https://www.rcsb.org/structure/4M8O) (Hogg et al. 2014)
- Roske & Yeeles, *Nat Struct Mol Biol* 31:1921-1931 (2024) — first human Pol e catalytic domain structures: open (9F6D), ajar (9F6E), and closed (9F6F) Finger conformations; four mismatch proofreading intermediates: post-insertion ([9F6I](https://www.rcsb.org/structure/9F6I)), arrest ([9F6J](https://www.rcsb.org/structure/9F6J)), frayed substrate ([9F6K](https://www.rcsb.org/structure/9F6K)), mismatch excision ([9F6L](https://www.rcsb.org/structure/9F6L))
- He, Wang, Yao, O'Donnell & Li, *Nat Commun* 15:7847 (2024) — independent human Pol e-PCNA cryo-EM structure (9B8S)
- Mur, Viana-Errasti et al., *Genome Med* 15:85 (2023) — gene-specific ACMG/AMP framework for POLE/POLD1 germline variant classification with likelihood-ratio-based Bayesian update
- Yuan, Georgescu, Schauer, O'Donnell & Li, *Nat Commun* 11:3156 (2020) — structure of the polymerase e holoenzyme and atomic model of the leading strand replisome
- Robinson, Coorens et al., *Nat Genet* 53:1434 (2021) — normal-tissue somatic mutation rates in POLE/POLD1 germline carriers (Sanger)
- Hogg, Osterman, Bylund, Ganai, Lundstrom, Sauer-Eriksson & Johansson, *Nat Struct Mol Biol* 21:49-55 (2014) — S. cerevisiae Pol2 catalytic domain (4M8O)
- Palles et al., *Nat Genet* 45:136-144 (2013) — germline POLE/POLD1 proofreading-domain mutations predispose to colorectal adenomas and carcinomas; defines PPAP and germline L424V
- Leon-Castillo et al., *J Pathol* 250:323-335 (2020) — POLE-score and 11 pathogenic somatic EDMs for endometrial carcinoma molecular classification
- Parkash et al., *Nat Commun* 10:373 (2019) — structural consequence of P286R (via yeast P301R); R blocks nascent-strand exo-site access
- ter Beek et al., *Nucleic Acids Res* 47:5712-5722 (2019) — structural evidence for [4Fe-4S] in the CysX motif of the Pol e catalytic core
- Jain et al., *J Mol Biol* 426:301-308 (2014) — iron-sulfur cluster in the polymerase domain of yeast Pol e
- Ganai, Bylund & Johansson, *Nucleic Acids Res* 43:932-942 (2015) — biochemical basis of pol<->exo shuttling
- Alexandrov et al., *Nature* 578:94-101 (2020) — mutational signatures SBS10a/10b/14/28 (POLE); note COSMIC v3.4 added SBS10c/SBS10d (POLD1-associated)
- Rayner et al., *Nat Rev Cancer* 16:71-81 (2016) — review of polymerase proofreading domain mutations in cancer

## License

[Apache 2.0](LICENSE) — see [SECURITY.md](SECURITY.md) for vulnerability reporting.
