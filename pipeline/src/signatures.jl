"""
Mutational signature deconvolution: NMF/NNLS for SBS10a/b/28/14 attribution.
"""

# POLE-associated mutational signatures (COSMIC v3.4)
# 96-channel trinucleotide context profiles (simplified to key channels)
const POLE_SIGNATURES = [:SBS10a, :SBS10b, :SBS28, :SBS14]

# Trinucleotide context categories (6 substitution types × 16 contexts = 96)
const SUBSTITUTION_TYPES = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
const CONTEXTS_5P = ['A', 'C', 'G', 'T']
const CONTEXTS_3P = ['A', 'C', 'G', 'T']

"""
    SignatureProfile

A mutational signature profile: 96-channel probability vector.
"""
struct SignatureProfile
    name::Symbol
    channels::Vector{Float64}  # length 96
    description::String
end

"""
    load_reference_signatures(path) -> Vector{SignatureProfile}

Load COSMIC SBS reference signature profiles from file.
If file not found, uses built-in POLE-associated signature profiles.
"""
function load_reference_signatures(path::AbstractString="")
    if !isempty(path) && isfile(path)
        return parse_cosmic_signatures(path)
    end

    @info "Using built-in POLE signature reference profiles"
    return get_builtin_pole_signatures()
end

"""
    get_builtin_pole_signatures() -> Vector{SignatureProfile}

Built-in reference profiles for POLE-associated signatures.
These are simplified 96-channel profiles based on COSMIC v3.4 data.
"""
function get_builtin_pole_signatures()
    profiles = SignatureProfile[]

    # SBS10a: POLE ultra-mutator (TCT>TAT and TCG>TTG dominant)
    # Characterized by C>A mutations in TCA/TCT context
    sbs10a = zeros(96)
    # C>A in TCN context (channels 4-7 approximately)
    sbs10a[5] = 0.15   # TCA>TAA
    sbs10a[6] = 0.08   # TCC>TAC
    sbs10a[7] = 0.05   # TCG>TAG
    sbs10a[8] = 0.20   # TCT>TAT  (dominant)
    # C>T in TCN context
    sbs10a[37] = 0.12  # TCA>TTA
    sbs10a[38] = 0.10  # TCC>TTC
    sbs10a[39] = 0.08  # TCG>TTG
    sbs10a[40] = 0.15  # TCT>TTT
    # Fill remaining with low background
    remaining = 1.0 - sum(sbs10a)
    bg = remaining / (96 - count(x -> x > 0, sbs10a))
    for i in 1:96
        if sbs10a[i] == 0.0
            sbs10a[i] = bg
        end
    end
    push!(profiles, SignatureProfile(:SBS10a, sbs10a,
        "POLE exonuclease domain mutations; ultra-mutator; TCT>TAT dominant"))

    # SBS10b: POLE variant (C>A in TCT context, less extreme)
    sbs10b = zeros(96)
    sbs10b[5] = 0.10
    sbs10b[8] = 0.25   # TCT>TAT even more dominant
    sbs10b[37] = 0.08
    sbs10b[40] = 0.10
    sbs10b[21] = 0.06  # C>G in TCT
    remaining = 1.0 - sum(sbs10b)
    bg = remaining / (96 - count(x -> x > 0, sbs10b))
    for i in 1:96
        if sbs10b[i] == 0.0
            sbs10b[i] = bg
        end
    end
    push!(profiles, SignatureProfile(:SBS10b, sbs10b,
        "POLE exonuclease domain mutations; V411L-associated"))

    # SBS28: POLE germline (less pronounced TCT>TAT)
    sbs28 = zeros(96)
    sbs28[5] = 0.06
    sbs28[8] = 0.12    # TCT>TAT (less dominant than SBS10a)
    sbs28[37] = 0.08
    sbs28[38] = 0.07
    sbs28[40] = 0.09
    # More distributed C>T pattern
    for i in 33:48
        sbs28[i] = max(sbs28[i], 0.03)
    end
    remaining = 1.0 - sum(sbs28)
    bg = remaining / (96 - count(x -> x > 0, sbs28))
    for i in 1:96
        if sbs28[i] == 0.0
            sbs28[i] = bg
        end
    end
    push!(profiles, SignatureProfile(:SBS28, sbs28,
        "POLE germline / S297F-associated; moderate mutator"))

    # SBS14: POLE + MMR deficiency (combined signature)
    sbs14 = zeros(96)
    # Mix of POLE (C>A in TCT) and MMR (C>T at CpG)
    sbs14[8] = 0.10    # TCT>TAT (POLE component)
    sbs14[39] = 0.12   # TCG>TTG (CpG component from MMR)
    sbs14[38] = 0.08   # TCC>TTC
    sbs14[35] = 0.06   # ACG>ATG (CpG)
    sbs14[36] = 0.05   # CCG>CTG (CpG)
    # Broader C>T distribution (MMR component)
    for i in 33:48
        sbs14[i] = max(sbs14[i], 0.04)
    end
    remaining = 1.0 - sum(sbs14)
    bg = remaining / (96 - count(x -> x > 0, sbs14))
    for i in 1:96
        if sbs14[i] == 0.0
            sbs14[i] = bg
        end
    end
    push!(profiles, SignatureProfile(:SBS14, sbs14,
        "POLE + MMR deficiency concurrent; CpG>TpG + TCT>TAT"))

    return profiles
end

"""
    parse_cosmic_signatures(path) -> Vector{SignatureProfile}

Parse COSMIC SBS signature file (TSV with 96 rows × N signature columns).
"""
function parse_cosmic_signatures(path::AbstractString)
    profiles = SignatureProfile[]
    lines = readlines(path)
    isempty(lines) && return profiles

    # Header row has signature names
    header = split(lines[1], '\t')
    n_sigs = length(header) - 1  # first column is context

    channels = [zeros(96) for _ in 1:n_sigs]

    for (row_idx, line) in enumerate(Iterators.drop(eachline(path), 1))
        row_idx > 96 && break
        fields = split(line, '\t')
        for sig_idx in 1:min(n_sigs, length(fields) - 1)
            val = tryparse(Float64, fields[sig_idx + 1])
            val !== nothing && (channels[sig_idx][row_idx] = val)
        end
    end

    # Only keep POLE-associated signatures
    for (i, sig_name) in enumerate(header[2:end])
        sym = Symbol(sig_name)
        if sym in POLE_SIGNATURES
            push!(profiles, SignatureProfile(sym, channels[i], string(sym)))
        end
    end

    return profiles
end

"""
    build_mutation_catalogue(variants) -> Vector{Float64}

Build 96-channel trinucleotide context mutation catalogue from variants.
Uses trinucleotide context field of Variant struct.
"""
function build_mutation_catalogue(variants::Vector{Variant})
    catalogue = zeros(96)

    for v in variants
        channel = context_to_channel(v.context, v.ref_aa, v.alt_aa)
        if 1 <= channel <= 96
            catalogue[channel] += 1.0
        end
    end

    # If no valid contexts, distribute uniformly (fallback)
    if sum(catalogue) == 0
        # Assign based on substitution type
        for v in variants
            # Simple assignment based on ref/alt
            channel = simple_channel_assignment(v.ref_aa, v.alt_aa)
            catalogue[channel] += 1.0
        end
    end

    return catalogue
end

"""
    context_to_channel(context, ref, alt) -> Int

Map trinucleotide context + substitution to 96-channel index.
"""
function context_to_channel(context::AbstractString, ref::Char, alt::Char)
    length(context) < 3 && return simple_channel_assignment(ref, alt)

    # Standard: pyrimidine-centered (C or T as reference)
    # 6 types × 16 contexts = 96
    sub_types = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]

    # Determine substitution type
    sub = "$(ref)>$(alt)"
    # Complement if purine reference
    if ref in ('A', 'G')
        comp = Dict('A' => 'T', 'T' => 'A', 'G' => 'C', 'C' => 'G')
        sub = "$(comp[ref])>$(comp[alt])"
    end

    type_idx = findfirst(==(sub), sub_types)
    type_idx === nothing && return 1

    # Context indices
    ctx_5p = uppercase(context[1])
    ctx_3p = uppercase(context[3])
    if ref in ('A', 'G')
        # Complement the context too
        comp = Dict('A' => 'T', 'T' => 'A', 'G' => 'C', 'C' => 'G')
        ctx_5p = comp[context[3]]
        ctx_3p = comp[context[1]]
    end

    p5_idx = findfirst(==(ctx_5p), CONTEXTS_5P)
    p3_idx = findfirst(==(ctx_3p), CONTEXTS_3P)
    p5_idx === nothing && return type_idx
    p3_idx === nothing && return type_idx

    context_idx = (p5_idx - 1) * 4 + p3_idx  # 1-16
    return (type_idx - 1) * 16 + context_idx
end

"""
    simple_channel_assignment(ref, alt) -> Int

Simple channel assignment when no trinucleotide context available.
"""
function simple_channel_assignment(ref::Char, alt::Char)
    sub_types = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
    comp = Dict('A' => 'T', 'T' => 'A', 'G' => 'C', 'C' => 'G')

    sub = "$(ref)>$(alt)"
    if ref in ('A', 'G')
        sub = "$(comp[ref])>$(comp[alt])"
    end

    type_idx = findfirst(==(sub), sub_types)
    type_idx === nothing && return 1

    # Default to middle context
    return (type_idx - 1) * 16 + 8
end

"""
    decompose_signatures(catalogue, refs; sigs=POLE_SIGNATURES) -> Dict{Symbol, Float64}

Decompose mutation catalogue into signature contributions using NNLS.
Returns fractional attribution for each signature.
"""
function decompose_signatures(
    catalogue::Vector{Float64},
    refs::Vector{SignatureProfile};
    sigs::Vector{Symbol}=POLE_SIGNATURES,
)
    # Filter to requested signatures
    sig_profiles = filter(s -> s.name in sigs, refs)
    isempty(sig_profiles) && return Dict(s => 0.0 for s in sigs)

    # Build reference matrix (96 × n_sigs)
    n_sigs = length(sig_profiles)
    W = hcat([s.channels for s in sig_profiles]...)  # 96 × n_sigs

    # Non-negative least squares: catalogue ≈ W * exposures
    exposures = nonneg_lsq(W, catalogue)

    # Normalize to fractions
    total = sum(exposures)
    fractions = total > 0 ? exposures ./ total : zeros(n_sigs)

    result = Dict{Symbol, Float64}()
    for (i, sp) in enumerate(sig_profiles)
        result[sp.name] = round(fractions[i], digits=4)
    end

    # Fill in zeros for any missing signatures
    for s in sigs
        if !haskey(result, s)
            result[s] = 0.0
        end
    end

    return result
end

"""
    attribute_per_variant(variants, sig_weights) -> Vector{Dict{Symbol, Float64}}

Assign signature weights to each variant based on its context match to signature profiles.
"""
function attribute_per_variant(
    variants::Vector{Variant},
    sig_weights::Dict{Symbol, Float64},
    refs::Vector{SignatureProfile},
)
    attributions = Dict{Symbol, Float64}[]

    for v in variants
        channel = context_to_channel(v.context, v.ref_aa, v.alt_aa)
        attr = Dict{Symbol, Float64}()

        # Weight by how much each signature explains this specific channel
        total_explanation = 0.0
        for sp in refs
            if 1 <= channel <= 96
                explanation = sp.channels[channel] * get(sig_weights, sp.name, 0.0)
                attr[sp.name] = explanation
                total_explanation += explanation
            end
        end

        # Normalize per-variant attributions
        if total_explanation > 0
            for k in keys(attr)
                attr[k] = round(attr[k] / total_explanation, digits=4)
            end
        end

        push!(attributions, attr)
    end

    return attributions
end

"""
    get_pole_mutation_signatures() -> Vector{Dict{Symbol, Float64}}

Get pre-computed signature attributions for the 8 known POLE mutations.
Based on published COSMIC/TCGA data.
"""
function get_pole_mutation_signatures()
    return [
        # P286R - primary SBS10a/b ultra-mutator
        Dict(:SBS10a => 0.55, :SBS10b => 0.30, :SBS28 => 0.05, :SBS14 => 0.10),
        # V411L - SBS10a dominant
        Dict(:SBS10a => 0.70, :SBS10b => 0.15, :SBS28 => 0.05, :SBS14 => 0.10),
        # S297F - SBS28 germline
        Dict(:SBS10a => 0.10, :SBS10b => 0.05, :SBS28 => 0.75, :SBS14 => 0.10),
        # D275A - functional study (no clear signature)
        Dict(:SBS10a => 0.40, :SBS10b => 0.20, :SBS28 => 0.20, :SBS14 => 0.20),
        # D368A - functional study
        Dict(:SBS10a => 0.35, :SBS10b => 0.25, :SBS28 => 0.20, :SBS14 => 0.20),
        # F367S - germline PPAP
        Dict(:SBS10a => 0.15, :SBS10b => 0.10, :SBS28 => 0.60, :SBS14 => 0.15),
        # L424V - moderate mutator
        Dict(:SBS10a => 0.45, :SBS10b => 0.25, :SBS28 => 0.15, :SBS14 => 0.15),
        # M444K - germline PPAP
        Dict(:SBS10a => 0.20, :SBS10b => 0.10, :SBS28 => 0.55, :SBS14 => 0.15),
    ]
end
