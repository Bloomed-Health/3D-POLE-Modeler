"""
DDG stability estimation for point mutants using statistical potential approach.
"""

# Amino acid hydrophobicity scale (Kyte-Doolittle, normalized to [-1, 1])
const AA_HYDROPHOBICITY = Dict{Char, Float64}(
    'A' =>  0.40,  'R' => -1.00,  'N' => -0.78,  'D' => -0.78,
    'C' =>  0.56,  'E' => -0.78,  'Q' => -0.78,  'G' =>  0.00,
    'H' => -0.67,  'I' =>  1.00,  'L' =>  0.84,  'K' => -0.87,
    'M' =>  0.42,  'F' =>  0.62,  'P' => -0.36,  'S' => -0.18,
    'T' => -0.05,  'W' => -0.20,  'Y' => -0.27,  'V' =>  0.91,
)

# Amino acid volume (Zamyatnin, normalized to [0, 1])
const AA_VOLUME = Dict{Char, Float64}(
    'A' => 0.25,  'R' => 0.78,  'N' => 0.43,  'D' => 0.40,
    'C' => 0.35,  'E' => 0.53,  'Q' => 0.55,  'G' => 0.00,
    'H' => 0.58,  'I' => 0.68,  'L' => 0.68,  'K' => 0.73,
    'M' => 0.63,  'F' => 0.75,  'P' => 0.38,  'S' => 0.28,
    'T' => 0.38,  'W' => 1.00,  'Y' => 0.83,  'V' => 0.55,
)

# Amino acid charge
const AA_FORMAL_CHARGE = Dict{Char, Float64}(
    'A' =>  0.0,  'R' =>  1.0,  'N' =>  0.0,  'D' => -1.0,
    'C' =>  0.0,  'E' => -1.0,  'Q' =>  0.0,  'G' =>  0.0,
    'H' =>  0.5,  'I' =>  0.0,  'L' =>  0.0,  'K' =>  1.0,
    'M' =>  0.0,  'F' =>  0.0,  'P' =>  0.0,  'S' =>  0.0,
    'T' =>  0.0,  'W' =>  0.0,  'Y' =>  0.0,  'V' =>  0.0,
)

# Secondary structure propensity differences (Chou-Fasman inspired)
const SS_PROPENSITY = Dict{Char, Dict{Symbol, Float64}}(
    'A' => Dict(:H => 1.42, :E => 0.83, :L => 0.66),
    'R' => Dict(:H => 0.98, :E => 0.93, :L => 1.01),
    'N' => Dict(:H => 0.67, :E => 0.89, :L => 1.56),
    'D' => Dict(:H => 1.01, :E => 0.54, :L => 1.46),
    'C' => Dict(:H => 0.70, :E => 1.19, :L => 1.19),
    'E' => Dict(:H => 1.51, :E => 0.37, :L => 0.74),
    'Q' => Dict(:H => 1.11, :E => 1.10, :L => 0.98),
    'G' => Dict(:H => 0.57, :E => 0.75, :L => 1.56),
    'H' => Dict(:H => 1.00, :E => 0.87, :L => 0.95),
    'I' => Dict(:H => 1.08, :E => 1.60, :L => 0.47),
    'L' => Dict(:H => 1.21, :E => 1.30, :L => 0.59),
    'K' => Dict(:H => 1.16, :E => 0.74, :L => 1.01),
    'M' => Dict(:H => 1.45, :E => 1.05, :L => 0.60),
    'F' => Dict(:H => 1.13, :E => 1.38, :L => 0.60),
    'P' => Dict(:H => 0.57, :E => 0.55, :L => 1.52),
    'S' => Dict(:H => 0.77, :E => 0.75, :L => 1.43),
    'T' => Dict(:H => 0.83, :E => 1.19, :L => 0.96),
    'W' => Dict(:H => 1.08, :E => 1.37, :L => 0.96),
    'Y' => Dict(:H => 0.69, :E => 1.47, :L => 1.14),
    'V' => Dict(:H => 1.06, :E => 1.70, :L => 0.50),
)

"""
    estimate_ddg(ref_aa, alt_aa, ss_type, burial, conservation) -> Float64

Estimate ΔΔG (kcal/mol) for a point mutation using statistical potential approach.

# Arguments
- `ref_aa::Char`: Wild-type amino acid (one-letter)
- `alt_aa::Char`: Mutant amino acid (one-letter)
- `ss_type::Symbol`: Secondary structure type (:H, :E, :L)
- `burial::Float64`: Residue burial score (0=surface, 1=buried)
- `conservation::Float64`: Conservation score (0=variable, 1=conserved)

# Returns
- DDG estimate in kcal/mol (positive = destabilizing, negative = stabilizing)

# Method
Combines:
1. Hydrophobicity mismatch × burial (buried hydrophobics stabilize)
2. Volume change × burial (steric clashes in core)
3. Charge change × (1-burial) (surface charge changes less costly)
4. SS propensity difference (helix breakers in helices, etc.)
5. Conservation penalty (mutations at conserved sites more costly)
"""
function estimate_ddg(
    ref_aa::Char, alt_aa::Char,
    ss_type::Symbol,
    burial::Float64,
    conservation::Float64,
)
    ref = uppercase(ref_aa)
    alt = uppercase(alt_aa)

    # 1. Hydrophobicity mismatch
    Δhydro = get(AA_HYDROPHOBICITY, alt, 0.0) - get(AA_HYDROPHOBICITY, ref, 0.0)
    # Buried: hydrophobic→polar is destabilizing; Surface: polar→hydrophobic less costly
    hydro_term = -Δhydro * burial * 1.5  # negative Δhydro in core = destabilizing

    # 2. Volume change (steric)
    Δvol = abs(get(AA_VOLUME, alt, 0.5) - get(AA_VOLUME, ref, 0.5))
    volume_term = Δvol * burial * 2.0  # larger changes in core more destabilizing

    # 3. Charge change
    Δcharge = abs(get(AA_FORMAL_CHARGE, alt, 0.0) - get(AA_FORMAL_CHARGE, ref, 0.0))
    charge_term = Δcharge * (1.0 - burial) * 0.5 + Δcharge * burial * 1.5

    # 4. SS propensity
    ref_prop = get(get(SS_PROPENSITY, ref, Dict()), ss_type, 1.0)
    alt_prop = get(get(SS_PROPENSITY, alt, Dict()), ss_type, 1.0)
    ss_term = max(0.0, ref_prop - alt_prop) * 1.0  # loss of propensity

    # 5. Conservation penalty
    cons_term = conservation * 1.5  # highly conserved → mutation more costly

    # Total DDG
    ddg = hydro_term + volume_term + charge_term + ss_term + cons_term

    # Scale to realistic range (most DDG values are -2 to +8 kcal/mol)
    ddg = clamp(ddg, -3.0, 10.0)

    return round(ddg, digits=2)
end

"""
    estimate_burial(residue_idx, residues; radius=10.0) -> Float64

Estimate residue burial from C-alpha contact density.
More contacts within radius → more buried.
"""
function estimate_burial(
    residue_idx::Int,
    residues::Vector{Residue3D};
    radius::Float64=10.0,
)
    target = residues[residue_idx]
    n_contacts = 0
    for (i, r) in enumerate(residues)
        i == residue_idx && continue
        d = sqrt((target.x - r.x)^2 + (target.y - r.y)^2 + (target.z - r.z)^2)
        if d <= radius
            n_contacts += 1
        end
    end

    # Normalize: typical buried residue has ~20-30 contacts, surface has ~5-10
    burial = clamp(n_contacts / 25.0, 0.0, 1.0)
    return burial
end

"""
    batch_ddg(mutations, residues, conservation_scores) -> Vector{Float64}

Compute DDG for all mutations, using structural context when available.
"""
function batch_ddg(mutations, residues::Vector{Residue3D}, conservation_scores::Vector{Float64})
    ddg_values = Float64[]

    for mut in mutations
        # Find residue in structure
        res_idx = findfirst(r -> r.resnum == mut.residue, residues)

        if res_idx !== nothing
            ss_type = residues[res_idx].ss_type
            burial = estimate_burial(res_idx, residues)
        else
            # Default values when residue not in structure
            ss_type = :L
            burial = 0.5
        end

        # Get conservation at this position
        cons = if mut.residue <= length(conservation_scores)
            conservation_scores[mut.residue]
        else
            0.5
        end

        ddg = estimate_ddg(mut.ref_aa, mut.alt_aa, ss_type, burial, cons)
        push!(ddg_values, ddg)
    end

    return ddg_values
end
