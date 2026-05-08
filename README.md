# 3D POLE Modeler

Interactive 3D structural viewer for the human DNA Polymerase epsilon (Pol e) catalytic subunit (p261, UniProt Q07864).

## Overview

This tool renders a scientifically accurate cartoon ribbon model of the POLE holoenzyme directly in the browser using Three.js (WebGL). It is designed for exploring the structural biology of Pol e, its catalytic mechanism, cancer-associated mutations, and the clinical consequences of rare germline variants.

## Features

- **Cartoon ribbon representation** — helices, beta-strands (with arrowheads), and loops rendered per secondary structure element across all six functional domains (NTD, EXO, Palm, Fingers, Thumb, CTD)
- **Two-metal-ion catalysis** — ball-and-stick active site with octahedral Mg2+ coordination, conserved aspartate residues (D640, D642, D860), and distance annotations
- **Incoming dNTP + Watson-Crick base pairing** — template nucleotide with hydrogen bond visualization
- **3'-to-5' exonuclease site** — two-metal proofreading active site with catalytic residues D275, E277, D368, D370
- **Zinc-finger motifs** — CysA/CysB Zn2+ tetrahedral coordination in the CTD
- **DNA double helix** — template and primer backbone tubes with base-pair rungs
- **Mutation hotspot markers** — P286R, V411L, S297F with interactive tooltips
- **Accessory subunits** — ghost representations of POLE2 (p59), POLE3 (p17), POLE4 (p12)
- **PCNA sliding clamp** — toroidal processivity factor with homotrimer symmetry markers
- **Beta-sheet hydrogen bonds** — inter-strand H-bond visualization
- **Exo-Pol shuttling path** — animated primer terminus transfer between active sites (~35 A)

### Color Modes

| Mode | Description |
|------|-------------|
| Domain | Default domain-based coloring |
| B-factor | Flexibility gradient (blue = rigid, red = flexible) |
| Conservation | Evolutionary conservation (blue = variable, purple = conserved) |
| Variant | COSMIC/ClinVar variant density |
| Charge | Electrostatic surface potential |
| Pathogenicity | Clinical significance spectrum (blue-green = benign, white = VUS, red = pathogenic) |

### Color Schemes

- **Default** — muted, publication-quality palette
- **PyMOL** — classic PyMOL-style bright colors
- **Accessible** — colorblind-friendly palette

## Mutation Viewer

A separate page (`mutation.html`) visualizes the ultra-rare **c.138del** germline variant (p.Leu46Phefs*8), showing:

- Side-by-side wild-type vs. truncated mutant comparison
- Domain retention bar chart (54 / 2,286 residues = 2.4% translated)
- Reading frame analysis with frameshift coloring

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
| `1`-`6` | Focus domain |
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

## References

- UniProt [Q07864](https://www.uniprot.org/uniprot/Q07864) (POLE_HUMAN)
- PDB: 4M8O, 6WJV, 8D32
- Lancey et al., *Nat Struct Mol Biol* 2020
- Hogg et al., *Nat Struct Mol Biol* 2014

## License

Proprietary - Bloomed Health
