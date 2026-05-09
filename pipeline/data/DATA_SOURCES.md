# Pipeline Data Sources

This directory holds input data for the Julia POLE pipeline. Files are NOT included
in the repository; download them from the sources below.

## Required Data Files

### PDB structures (`raw/`)
Downloaded automatically by the pipeline from RCSB PDB.
- 9F6D (Roske & Yeeles 2024 — open Finger)
- 9F6E (ajar)
- 9F6F (closed)
- 9B8S (He et al. 2024)
- 4M8O (Hogg et al. 2014 — yeast Pol2)

### ClinVar variants (`variants/clinvar_pole.tsv`)
1. Go to https://www.ncbi.nlm.nih.gov/clinvar/
2. Search: `POLE[gene] AND "single nucleotide variant"[vartype]`
3. Download: Send to > File > Tab-delimited
4. Expected columns: Gene, Name (HGVS), Clinical significance, Review status

### COSMIC somatic variants (`variants/cosmic_pole.tsv`)
1. Register at https://cancer.sanger.ac.uk/cosmic/register
2. Download COSMIC Mutation Data (requires academic license)
3. Filter for POLE gene
4. Expected columns: Gene, AA Mutation, Count, Primary site

### AlphaMissense predictions (`variants/alphamissense_pole.tsv`)
1. Download from https://zenodo.org/records/8208688 (Cheng et al. 2023)
2. Filter for UniProt Q07864 (POLE_HUMAN)
3. Expected columns: uniprot_id, position, ref_aa, alt_aa, am_pathogenicity, am_class

### POLE ortholog MSA (`msa/pole_orthologs.fasta`)
For per-residue conservation (ConSurf-style):
1. Go to https://www.uniprot.org/uniref/UniRef50_Q07864
2. Download aligned sequences as FASTA
3. Or run ConSurf: https://consurf.tau.ac.il/ with UniProt Q07864
4. Save per-residue scores as FASTA alignment

### COSMIC SBS signatures (`signatures/COSMIC_v3.4_SBS_GRCh38.txt`)
1. Download from https://cancer.sanger.ac.uk/signatures/documents/2123/COSMIC_v3.4_SBS_GRCh38.txt
2. Official 96-channel mutational signature profiles
3. Reference: Alexandrov et al., Nature 578:94-101 (2020)

## Fallback Behavior

When data files are absent, the pipeline falls back to built-in values:
- **ClinVar/COSMIC**: Uses ~15 known pathogenic/likely pathogenic POLE variants
- **AlphaMissense**: Skipped (no fallback)
- **Conservation**: Generates synthetic per-domain estimates (NOT real per-residue data)
- **Signatures**: Uses hand-curated approximations of SBS10a/b/28/14 profiles
- **PDB**: Downloads from RCSB automatically

Built-in fallback values are clearly labeled in the JSON output with
`"source": "built_in"` vs `"source": "file"` provenance tracking.
