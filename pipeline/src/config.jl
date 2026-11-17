"""
Pipeline configuration: PDB IDs, domain ranges, output paths.
"""

# Domain boundary residue ranges (UniProt Q07864 / POLE catalytic subunit)
const DOMAIN_RANGES = OrderedDict(
    :ntd      => (1, 222),
    :exo      => (223, 517),
    :palm     => (518, 870),
    :pdomain  => (640, 760),
    :fingers  => (870, 1020),
    :thumb    => (1020, 1189),
    :inactpol => (1190, 1900),
    :ctd      => (1900, 2286),
)

const DOMAIN_NAMES = OrderedDict(
    :ntd      => "N-terminal domain",
    :exo      => "3'→5' Exonuclease",
    :palm     => "Palm subdomain",
    :pdomain  => "Processivity (P) domain",
    :fingers  => "Fingers subdomain",
    :thumb    => "Thumb subdomain",
    :inactpol => "Inactive Pol module",
    :ctd      => "C-terminal domain",
)

const DOMAIN_ABBREVS = OrderedDict(
    :ntd      => "NTD",
    :exo      => "EXO",
    :palm     => "PALM",
    :pdomain  => "P-DOM",
    :fingers  => "FNG",
    :thumb    => "THM",
    :inactpol => "INACT",
    :ctd      => "CTD",
)

const TOTAL_RESIDUES = 2286

# PDB entries for POLE structures
const PDB_IDS = ["4M8O", "9F6D", "9F6E", "9F6F", "9B8S"]
const PRIMARY_PDB = "9F6D"  # Human Pol epsilon (Roske & Yeeles 2024)

# Known POLE pathogenic mutations
# Clinical and somatic POLE mutations (viewer markers)
const POLE_MUTATIONS = [
    (id="P286R",  residue=286,  ref_aa='P', alt_aa='R', domain=:exo,
     label="Pro286→Arg",
     detail="Somatic ultra-mutator / adjacent to ExoI / CRC & endometrial / SBS10a/b"),
    (id="V411L",  residue=411,  ref_aa='V', alt_aa='L', domain=:exo,
     label="Val411→Leu",
     detail="Somatic/germline proofreading-deficient / TMB >100 mut/Mb / SBS10a; Barbari 2018 (G3 yeast); human data Hodel 2020; REVEL 0.457"),
    (id="S297F",  residue=297,  ref_aa='S', alt_aa='F', domain=:exo,
     label="Ser297→Phe",
     detail="Somatic hotspot / outside Exo motifs / endometrial & CRC / REVEL 0.799"),
    (id="L424V",  residue=424,  ref_aa='L', alt_aa='V', domain=:exo,
     label="Leu424→Val",
     detail="Most recurrent germline PPAP (24 families; Palles 2013) / within Exo IV / Note: somatic L424I is distinct"),
    (id="D287E",  residue=287,  ref_aa='D', alt_aa='E', domain=:exo,
     label="Asp287→Glu",
     detail="Benign population variant (gnomAD MAF ~0.08%, 216/268316 alleles; REVEL 0.286)"),
    (id="P436R",  residue=436,  ref_aa='P', alt_aa='R', domain=:exo,
     label="Pro436→Arg",
     detail="Somatic/germline Exo V / endometrial / SBS10a; germline CMMRD-like cases reported"),
    (id="M444K",  residue=444,  ref_aa='M', alt_aa='K', domain=:exo,
     label="Met444→Lys",
     detail="Germline PPAP / flanks Exo V / posterior 0.88 < 0.90 threshold"),
    (id="S459F",  residue=459,  ref_aa='S', alt_aa='F', domain=:exo,
     label="Ser459→Phe",
     detail="Somatic hotspot / endometrial / Leon-Castillo 2020 EDM set"),
    (id="F367S",  residue=367,  ref_aa='F', alt_aa='S', domain=:exo,
     label="Phe367→Ser",
     detail="Somatic (Leon-Castillo 2020) / in Exo II / c.1100T>C"),
    (id="A456P",  residue=456,  ref_aa='A', alt_aa='P', domain=:exo,
     label="Ala456→Pro",
     detail="Somatic EDM / inside Exo III / 19 tumours / Leon-Castillo 2020"),
    (id="M295R",  residue=295,  ref_aa='M', alt_aa='R', domain=:exo,
     label="Met295→Arg",
     detail="Somatic EDM / Leon-Castillo 2020 / endometrial"),
    (id="L424I",  residue=424,  ref_aa='L', alt_aa='I', domain=:exo,
     label="Leu424→Ile",
     detail="Somatic EDM / Leon-Castillo 2020 / distinct from germline L424V"),
    (id="D368Y",  residue=368,  ref_aa='D', alt_aa='Y', domain=:exo,
     label="Asp368→Tyr",
     detail="Somatic EDM / ExoII catalytic D368 / Leon-Castillo 2020"),
]

# Frameshift variant data (c.138del)
const FRAMESHIFT_VARIANT = (
    hgvs_c = "c.138del",
    hgvs_p = "p.Leu46Phefs*8",
    truncation_site = 53,  # stop codon at position 53; 52 residues translated
    wt_codons = [
        (codon="GAG", aa="Glu", pos=43),
        (codon="AGC", aa="Ser", pos=44),
        (codon="TCC", aa="Ser", pos=45),
        (codon="CTG", aa="Leu", pos=46),
        (codon="GAG", aa="Glu", pos=47),
        (codon="TTG", aa="Leu", pos=48),
        (codon="AGG", aa="Arg", pos=49),
        (codon="GCA", aa="Ala", pos=50),
        (codon="GCT", aa="Ala", pos=51),
        (codon="TTA", aa="Leu", pos=52),
        (codon="CCT", aa="Pro", pos=53),
    ],
    mut_codons = [
        (codon="GAG", aa="Glu", pos=43, status="normal"),
        (codon="AGC", aa="Ser", pos=44, status="normal"),
        (codon="TCC", aa="Ser", pos=45, status="normal"),
        (codon="TGA", aa="Phe", pos=46, status="frameshift"),
        (codon="GGA", aa="Gly", pos=47, status="frameshift"),
        (codon="GTT", aa="Val", pos=48, status="frameshift"),
        (codon="GAG", aa="Glu", pos=49, status="frameshift"),
        (codon="GGC", aa="Gly", pos=50, status="frameshift"),
        (codon="AGC", aa="Ser", pos=51, status="frameshift"),
        (codon="TTT", aa="Phe", pos=52, status="frameshift"),
        (codon="TAG", aa="*",   pos=53, status="stop"),
    ],
)

# Catalytic residues
const CATALYTIC_RESIDUES = Dict(
    :exo_site => [275, 277, 368, 462],      # D275, E277 (ExoI), D368 (ExoII, DEDDy catalytic set), D462 (ExoIII, verified from PDB 9F6L)
    :pol_motif_a => [626],                   # D626 (He et al. 2024, 9B8S)
    :pol_motif_c => [860, 862],              # D860, D862 (He et al. 2024, 9B8S)
)

# Metal sites
const METAL_TYPES = ["MG", "ZN", "FE", "S"]
const METAL_COORDINATION_CUTOFF = 2.5  # Angstroms

# AlphaFold
const ALPHAFOLD_UNIPROT = "Q07864"
const ALPHAFOLD_MODEL_VERSION = 6

"""
    PipelineConfig

Configuration for a pipeline run.
"""
struct PipelineConfig
    pdb_ids::Vector{String}
    primary_pdb::String
    primary_chain::String
    domain_ranges::OrderedDict{Symbol, Tuple{Int,Int}}
    total_residues::Int
    output_dir::String
    data_dir::String
    cache_dir::String
    conservation_window::Int
    density_window::Int
    charge_window::Int
end

function PipelineConfig(;
    pdb_ids = PDB_IDS,
    primary_pdb = PRIMARY_PDB,
    primary_chain = "A",
    domain_ranges = DOMAIN_RANGES,
    total_residues = TOTAL_RESIDUES,
    output_dir = joinpath(@__DIR__, "..", "..", "data"),
    data_dir = joinpath(@__DIR__, "..", "data"),
    cache_dir = joinpath(@__DIR__, "..", "data", "raw"),
    conservation_window = 5,
    density_window = 10,
    charge_window = 5,
)
    PipelineConfig(
        pdb_ids, primary_pdb, primary_chain,
        domain_ranges, total_residues,
        output_dir, data_dir, cache_dir,
        conservation_window, density_window, charge_window,
    )
end
