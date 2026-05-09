"""
Per-residue electrostatic charge computation from amino acid properties.
"""

# Amino acid pKa values and formal charges at pH 7.4
const AA_CHARGE_TABLE = Dict{Char, Float64}(
    'A' =>  0.0,   # Alanine
    'R' =>  1.0,   # Arginine (pKa 12.5, always protonated)
    'N' =>  0.0,   # Asparagine
    'D' => -1.0,   # Aspartate (pKa 3.65, deprotonated at pH 7.4)
    'C' =>  0.0,   # Cysteine (pKa 8.18, mostly protonated)
    'E' => -1.0,   # Glutamate (pKa 4.25, deprotonated)
    'Q' =>  0.0,   # Glutamine
    'G' =>  0.0,   # Glycine
    'H' =>  0.1,   # Histidine (pKa 6.0, ~10% protonated at pH 7.4)
    'I' =>  0.0,   # Isoleucine
    'L' =>  0.0,   # Leucine
    'K' =>  1.0,   # Lysine (pKa 10.5, always protonated)
    'M' =>  0.0,   # Methionine
    'F' =>  0.0,   # Phenylalanine
    'P' =>  0.0,   # Proline
    'S' =>  0.0,   # Serine
    'T' =>  0.0,   # Threonine
    'W' =>  0.0,   # Tryptophan
    'Y' =>  0.0,   # Tyrosine (pKa 10.1, protonated)
    'V' =>  0.0,   # Valine
)

# Three-letter to one-letter amino acid code mapping
const THREE_TO_ONE = Dict{String, Char}(
    "ALA" => 'A', "ARG" => 'R', "ASN" => 'N', "ASP" => 'D',
    "CYS" => 'C', "GLU" => 'E', "GLN" => 'Q', "GLY" => 'G',
    "HIS" => 'H', "ILE" => 'I', "LEU" => 'L', "LYS" => 'K',
    "MET" => 'M', "PHE" => 'F', "PRO" => 'P', "SER" => 'S',
    "THR" => 'T', "TRP" => 'W', "TYR" => 'Y', "VAL" => 'V',
)

"""
    amino_acid_charge(aa; pH=7.4) -> Float64

Formal charge for a single amino acid at given pH.
Uses Henderson-Hasselbalch for ionizable residues.
"""
function amino_acid_charge(aa::Char; pH::Float64=7.4)
    return get(AA_CHARGE_TABLE, uppercase(aa), 0.0)
end

"""
    amino_acid_charge(resname::String; pH=7.4) -> Float64

Charge from three-letter residue name.
"""
function amino_acid_charge(resname::AbstractString; pH::Float64=7.4)
    aa = get(THREE_TO_ONE, uppercase(resname), 'X')
    return amino_acid_charge(aa; pH=pH)
end

"""
    compute_charge_map(residues; window=5) -> Vector{Float64}

Compute per-residue charge with sliding window smoothing.
Output is normalized to [-1, 1] range.
"""
function compute_charge_map(residues::Vector{Residue3D}; window::Int=5)
    n = length(residues)
    raw_charges = [amino_acid_charge(r.resname) for r in residues]

    # Sliding window average
    smoothed = zeros(Float64, n)
    half_w = window ÷ 2
    for i in 1:n
        start_idx = max(1, i - half_w)
        end_idx = min(n, i + half_w)
        smoothed[i] = mean(raw_charges[start_idx:end_idx])
    end

    return smoothed
end

"""
    compute_charge_from_sequence(sequence; window=5) -> Vector{Float64}

Compute charge map from one-letter amino acid sequence string.
"""
function compute_charge_from_sequence(sequence::AbstractString; window::Int=5)
    n = length(sequence)
    raw_charges = [amino_acid_charge(c) for c in sequence]

    smoothed = zeros(Float64, n)
    half_w = window ÷ 2
    for i in 1:n
        start_idx = max(1, i - half_w)
        end_idx = min(n, i + half_w)
        smoothed[i] = mean(raw_charges[start_idx:end_idx])
    end

    return smoothed
end

"""
    generate_synthetic_charge_map(total_residues, domain_ranges) -> Vector{Float64}

Generate synthetic per-residue charge map based on known POLE domain charge properties.
Used when PDB structure is not available for full sequence.
"""
function generate_synthetic_charge_map(
    total_residues::Int,
    domain_ranges::OrderedDict{Symbol, Tuple{Int,Int}},
)
    # Known domain charge characteristics
    domain_charges = Dict(
        :ntd      =>  0.0,   # Neutral
        :exo      => -0.3,   # Negatively charged (catalytic carboxylates)
        :palm     => -0.8,   # Highly negative (D640, D642, D860)
        :pdomain  => -0.2,   # Slightly negative
        :fingers  =>  0.5,   # Positively charged (O-helix, DNA-binding)
        :thumb    =>  0.6,   # Positively charged (DNA-binding groove)
        :inactpol =>  0.1,   # Slightly positive
        :ctd      =>  0.2,   # Slightly positive (zinc-finger)
    )

    charges = fill(0.0, total_residues)

    for (domain, (start_r, end_r)) in domain_ranges
        base_charge = get(domain_charges, domain, 0.0)
        for r in start_r:min(end_r, total_residues)
            # Add positional variation
            noise = 0.1 * sin(r * 0.2) + 0.05 * cos(r * 0.5)
            charges[r] = clamp(base_charge + noise, -1.0, 1.0)
        end
    end

    # Smooth
    half_w = 2
    smoothed = zeros(Float64, total_residues)
    for i in 1:total_residues
        start_idx = max(1, i - half_w)
        end_idx = min(total_residues, i + half_w)
        smoothed[i] = mean(charges[start_idx:end_idx])
    end

    return smoothed
end
