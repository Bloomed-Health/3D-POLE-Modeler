"""
Pipeline configuration: PDB IDs, domain ranges, output paths.
"""

# Domain boundary residue ranges (UniProt Q07864 / POLE catalytic subunit)
const DOMAIN_RANGES = OrderedDict(
    :ntd      => (1, 280),
    :exo      => (268, 471),
    :palm     => (600, 870),
    :pdomain  => (700, 760),
    :fingers  => (870, 1020),
    :thumb    => (1020, 1190),
    :inactpol => (1198, 1900),
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
const PDB_IDS = ["4M8O", "6WJV", "8D32", "9F6D", "9F6E", "9F6F", "9B8S"]
const PRIMARY_PDB = "4M8O"  # Primary structure for initial pipeline

# Known POLE pathogenic mutations
const POLE_MUTATIONS = [
    (id="P286R",  residue=286,  ref_aa='P', alt_aa='R', domain=:exo,
     label="Pro286→Arg",
     detail="Somatic ultra-mutator / ExoII motif / CRC & endometrial / SBS10a/b"),
    (id="V411L",  residue=411,  ref_aa='V', alt_aa='L', domain=:exo,
     label="Val411→Leu",
     detail="Somatic proofreading-deficient / TMB >100 mut/Mb / SBS10a"),
    (id="S297F",  residue=297,  ref_aa='S', alt_aa='F', domain=:exo,
     label="Ser297→Phe",
     detail="Germline ExoI / PPAP predisposition / SBS28"),
    (id="D275A",  residue=275,  ref_aa='D', alt_aa='A', domain=:exo,
     label="Asp275→Ala",
     detail="ExoI catalytic / abolishes proofreading"),
    (id="D368A",  residue=368,  ref_aa='D', alt_aa='A', domain=:exo,
     label="Asp368→Ala",
     detail="ExoIII catalytic / Mg²⁺ coordination"),
    (id="F367S",  residue=367,  ref_aa='F', alt_aa='S', domain=:exo,
     label="Phe367→Ser",
     detail="Germline PPAP / adjacent to catalytic triad"),
    (id="L424V",  residue=424,  ref_aa='L', alt_aa='V', domain=:exo,
     label="Leu424→Val",
     detail="Somatic / moderate mutator phenotype"),
    (id="M444K",  residue=444,  ref_aa='M', alt_aa='K', domain=:exo,
     label="Met444→Lys",
     detail="Germline PPAP / exonuclease domain"),
]

# Frameshift variant data (c.138del)
const FRAMESHIFT_VARIANT = (
    hgvs_c = "c.138del",
    hgvs_p = "p.Glu47Argfs*8",
    truncation_site = 54,
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
        (codon="TAG", aa="*",   pos=54),
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
        (codon="ACC", aa="Thr", pos=53, status="frameshift"),
        (codon="TAG", aa="*",   pos=54, status="stop"),
    ],
)

# Catalytic residues
const CATALYTIC_RESIDUES = Dict(
    :exo_site => [275, 277, 368, 370],      # D275, E277, D368, D370
    :pol_motif_a => [640, 642],              # D640, D642
    :pol_motif_c => [860],                   # D860
)

# Metal sites
const METAL_TYPES = ["MG", "ZN", "FE", "S"]
const METAL_COORDINATION_CUTOFF = 2.5  # Angstroms

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
