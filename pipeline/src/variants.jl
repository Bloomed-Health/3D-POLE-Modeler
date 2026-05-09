"""
Variant database processing: ClinVar, COSMIC, AlphaMissense loading and scoring.
"""

"""
    Variant

A protein variant with source annotations.
"""
struct Variant
    residue::Int
    ref_aa::Char
    alt_aa::Char
    source::Symbol          # :clinvar, :cosmic, :alphamissense
    classification::String  # pathogenic, likely_pathogenic, VUS, etc.
    score::Float64          # pathogenicity score (0-1)
    frequency::Float64      # population/somatic frequency
    context::String         # trinucleotide context (for signatures)
end

"""
    load_clinvar_variants(path, gene="POLE") -> Vector{Variant}

Load ClinVar variants from TSV export.
Expected columns: gene, hgvs_p, classification, review_status
"""
function load_clinvar_variants(path::AbstractString; gene::String="POLE")
    variants = Variant[]

    if !isfile(path)
        @warn "ClinVar file not found: $path, using built-in POLE variants"
        return get_known_clinvar_variants()
    end

    for line in Iterators.drop(eachline(path), 1)  # skip header
        fields = split(line, '\t')
        length(fields) < 4 && continue
        fields[1] != gene && continue

        hgvs = fields[2]
        classif = fields[3]

        # Parse p.X###Y format
        m = match(r"p\.([A-Z])(\d+)([A-Z])", hgvs)
        m === nothing && continue

        ref = first(m.captures[1])
        pos = parse(Int, m.captures[2])
        alt = first(m.captures[3])

        score = classification_to_score(classif)
        push!(variants, Variant(pos, ref, alt, :clinvar, classif, score, 0.0, ""))
    end

    @info "Loaded $(length(variants)) ClinVar variants"
    return variants
end

"""
    load_cosmic_variants(path, gene="POLE") -> Vector{Variant}

Load COSMIC somatic variants from TSV.
Expected columns: gene, aa_mutation, count, primary_site
"""
function load_cosmic_variants(path::AbstractString; gene::String="POLE")
    variants = Variant[]

    if !isfile(path)
        @warn "COSMIC file not found: $path, using built-in somatic variants"
        return get_known_cosmic_variants()
    end

    for line in Iterators.drop(eachline(path), 1)
        fields = split(line, '\t')
        length(fields) < 4 && continue
        fields[1] != gene && continue

        aa_mut = fields[2]
        count = tryparse(Int, fields[3])
        count === nothing && (count = 1)

        m = match(r"p\.([A-Z])(\d+)([A-Z])", aa_mut)
        m === nothing && continue

        ref = first(m.captures[1])
        pos = parse(Int, m.captures[2])
        alt = first(m.captures[3])

        # Frequency proxy from recurrence count
        freq = min(count / 100.0, 1.0)
        score = min(0.5 + freq * 0.5, 1.0)  # recurrent = more likely pathogenic

        push!(variants, Variant(pos, ref, alt, :cosmic, "somatic", score, freq, ""))
    end

    @info "Loaded $(length(variants)) COSMIC variants"
    return variants
end

"""
    load_alphamissense_scores(path, gene="POLE") -> Vector{Variant}

Load AlphaMissense pathogenicity predictions.
Expected columns: uniprot_id, position, ref_aa, alt_aa, am_pathogenicity, am_class
"""
function load_alphamissense_scores(path::AbstractString; gene::String="POLE")
    variants = Variant[]

    if !isfile(path)
        @warn "AlphaMissense file not found: $path"
        return variants
    end

    for line in Iterators.drop(eachline(path), 1)
        fields = split(line, '\t')
        length(fields) < 6 && continue

        pos = tryparse(Int, fields[2])
        pos === nothing && continue
        ref = first(fields[3])
        alt = first(fields[4])
        score = tryparse(Float64, fields[5])
        score === nothing && continue
        classif = fields[6]

        push!(variants, Variant(pos, ref, alt, :alphamissense, classif, score, 0.0, ""))
    end

    @info "Loaded $(length(variants)) AlphaMissense scores"
    return variants
end

"""
    get_known_clinvar_variants() -> Vector{Variant}

Built-in known pathogenic POLE ClinVar variants when file not available.
"""
function get_known_clinvar_variants()
    return [
        Variant(286, 'P', 'R', :clinvar, "Pathogenic", 0.99, 0.0, "TCT>TGT"),
        Variant(411, 'V', 'L', :clinvar, "Pathogenic", 0.95, 0.0, "GTG>CTG"),
        Variant(297, 'S', 'F', :clinvar, "Pathogenic", 0.92, 0.0, "TCC>TTC"),
        Variant(275, 'D', 'A', :clinvar, "Pathogenic", 0.90, 0.0, "GAT>GCT"),
        Variant(368, 'D', 'A', :clinvar, "Pathogenic", 0.90, 0.0, "GAT>GCT"),
        Variant(367, 'F', 'S', :clinvar, "Likely_pathogenic", 0.85, 0.0, "TTC>TCC"),
        Variant(424, 'L', 'V', :clinvar, "Likely_pathogenic", 0.80, 0.0, "CTG>GTG"),
        Variant(444, 'M', 'K', :clinvar, "Likely_pathogenic", 0.78, 0.0, "ATG>AAG"),
        # Additional known variants
        Variant(396, 'N', 'D', :clinvar, "VUS", 0.50, 0.0, ""),
        Variant(459, 'P', 'L', :clinvar, "VUS", 0.45, 0.0, ""),
        Variant(640, 'D', 'N', :clinvar, "Pathogenic", 0.95, 0.0, "GAC>AAC"),
        Variant(642, 'D', 'N', :clinvar, "Pathogenic", 0.93, 0.0, "GAC>AAC"),
        Variant(860, 'D', 'N', :clinvar, "Pathogenic", 0.93, 0.0, "GAC>AAC"),
    ]
end

"""
    get_known_cosmic_variants() -> Vector{Variant}

Built-in known COSMIC somatic POLE variants.
"""
function get_known_cosmic_variants()
    return [
        Variant(286, 'P', 'R', :cosmic, "somatic", 0.99, 0.85, "TCT>TGT"),
        Variant(411, 'V', 'L', :cosmic, "somatic", 0.95, 0.45, "GTG>CTG"),
        Variant(297, 'S', 'F', :cosmic, "somatic", 0.90, 0.20, "TCC>TTC"),
        Variant(286, 'P', 'H', :cosmic, "somatic", 0.88, 0.10, "CCC>CAC"),
        Variant(286, 'P', 'L', :cosmic, "somatic", 0.85, 0.05, "CCT>CTT"),
        Variant(459, 'P', 'L', :cosmic, "somatic", 0.60, 0.08, "CCT>CTT"),
        Variant(367, 'F', 'S', :cosmic, "somatic", 0.82, 0.12, "TTC>TCC"),
        Variant(424, 'L', 'V', :cosmic, "somatic", 0.75, 0.15, "CTG>GTG"),
    ]
end

"""
    classification_to_score(classif) -> Float64

Convert ClinVar classification string to numeric score.
"""
function classification_to_score(classif::AbstractString)
    c = lowercase(classif)
    if contains(c, "pathogenic") && !contains(c, "likely")
        return 0.95
    elseif contains(c, "likely_pathogenic") || contains(c, "likely pathogenic")
        return 0.80
    elseif contains(c, "uncertain") || contains(c, "vus")
        return 0.50
    elseif contains(c, "likely_benign") || contains(c, "likely benign")
        return 0.20
    elseif contains(c, "benign") && !contains(c, "likely")
        return 0.05
    else
        return 0.50  # default to uncertain
    end
end

"""
    compute_variant_density(variants, total_residues; window=10) -> Vector{Float64}

Compute per-residue variant density using sliding window, normalized to [0, 1].
"""
function compute_variant_density(
    variants::Vector{Variant},
    total_residues::Int;
    window::Int=10,
)
    # Count variants per position
    counts = zeros(Float64, total_residues)
    for v in variants
        if 1 <= v.residue <= total_residues
            counts[v.residue] += 1.0
        end
    end

    # Sliding window sum
    density = zeros(Float64, total_residues)
    half_w = window ÷ 2
    for i in 1:total_residues
        start_idx = max(1, i - half_w)
        end_idx = min(total_residues, i + half_w)
        density[i] = sum(counts[start_idx:end_idx])
    end

    # Normalize to [0, 1]
    max_density = maximum(density)
    if max_density > 0
        density ./= max_density
    end

    return density
end

"""
    compute_pathogenicity_map(clinvar, cosmic, alphamissense, total_residues) -> Vector{Float64}

Combine variant sources into per-residue pathogenicity score.
Weights: ClinVar (0.5), COSMIC recurrence (0.3), AlphaMissense (0.2).
"""
function compute_pathogenicity_map(
    clinvar::Vector{Variant},
    cosmic::Vector{Variant},
    alphamissense::Vector{Variant},
    total_residues::Int,
)
    scores = fill(0.05, total_residues)  # baseline

    # ClinVar - highest weight
    for v in clinvar
        if 1 <= v.residue <= total_residues
            scores[v.residue] = max(scores[v.residue], v.score * 0.5 + scores[v.residue] * 0.5)
        end
    end

    # COSMIC recurrence
    for v in cosmic
        if 1 <= v.residue <= total_residues
            scores[v.residue] = max(scores[v.residue], v.score * 0.3 + scores[v.residue] * 0.7)
        end
    end

    # AlphaMissense
    for v in alphamissense
        if 1 <= v.residue <= total_residues
            scores[v.residue] = max(scores[v.residue], v.score * 0.2 + scores[v.residue] * 0.8)
        end
    end

    return clamp.(scores, 0.0, 1.0)
end
