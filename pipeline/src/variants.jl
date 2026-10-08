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

Placeholder: returns an empty variant list. Real ClinVar data should be fetched
from the ClinVar API or a downloaded VCF. The previous implementation contained
invented pathogenicity scores that could be mistaken for real ClinVar data.

To populate: use the ClinVar E-utilities API or download the ClinVar VCF for POLE
(Gene ID 5426) and parse with an appropriate VCF library.
"""
function get_known_clinvar_variants()
    @warn "get_known_clinvar_variants() returns empty list — real ClinVar data not yet integrated"
    return Variant[]
end

"""
    get_known_cosmic_variants() -> Vector{Variant}

Placeholder: returns an empty variant list. Real COSMIC data requires a licence
and should be loaded from a downloaded COSMIC VCF or TSV file.

To populate: download the COSMIC coding mutations file for POLE and parse.
"""
function get_known_cosmic_variants()
    @warn "get_known_cosmic_variants() returns empty list — real COSMIC data not yet integrated"
    return Variant[]
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
    density = sliding_window(counts, window; op=sum)

    # Normalize to [0, 1]
    max_density = maximum(density)
    if max_density > 0
        density ./= max_density
    end

    return density
end

"""
Pathogenicity source weights for combining ClinVar, COSMIC, and AlphaMissense scores.

Rationale:
- ClinVar (0.5): Expert-reviewed clinical variant classifications with structured
  evidence codes. Highest weight because assessments are curated by clinical labs
  and reviewed against ACMG/AMP guidelines.
- COSMIC (0.3): Somatic mutation recurrence across cancer cohorts. Recurrent somatic
  mutations in a specific gene/domain are strong evidence for functional impact,
  but lack germline clinical interpretation.
- AlphaMissense (0.2): Computational pathogenicity predictor based on protein language
  models. Genome-wide coverage but lower weight because purely computational with
  no clinical review.

Note: These weights are heuristic and not derived from a formal calibration study.
For clinical use, individual source scores should be evaluated independently.
"""
const PATHOGENICITY_WEIGHTS = (clinvar=0.5, cosmic=0.3, alphamissense=0.2)

"""
    compute_pathogenicity_map(clinvar, cosmic, alphamissense, total_residues;
                              weights=PATHOGENICITY_WEIGHTS) -> Vector{Float64}

Combine variant sources into per-residue pathogenicity score.
Weights: ClinVar (0.5), COSMIC recurrence (0.3), AlphaMissense (0.2).

The weights are heuristic: ClinVar = expert-reviewed clinical classifications,
COSMIC = somatic recurrence across cancer cohorts, AlphaMissense = computational
protein language model predictions. See `PATHOGENICITY_WEIGHTS` for rationale.
"""
function compute_pathogenicity_map(
    clinvar::Vector{Variant},
    cosmic::Vector{Variant},
    alphamissense::Vector{Variant},
    total_residues::Int;
    weights::NamedTuple=PATHOGENICITY_WEIGHTS,
)
    scores = fill(0.05, total_residues)  # baseline

    w_cv = weights.clinvar
    w_co = weights.cosmic
    w_am = weights.alphamissense

    # ClinVar - highest weight (expert-reviewed clinical classifications)
    for v in clinvar
        if 1 <= v.residue <= total_residues
            scores[v.residue] = max(scores[v.residue], v.score * w_cv + scores[v.residue] * (1 - w_cv))
        end
    end

    # COSMIC recurrence (somatic mutation frequency across cancer cohorts)
    for v in cosmic
        if 1 <= v.residue <= total_residues
            scores[v.residue] = max(scores[v.residue], v.score * w_co + scores[v.residue] * (1 - w_co))
        end
    end

    # AlphaMissense (computational pathogenicity predictor)
    for v in alphamissense
        if 1 <= v.residue <= total_residues
            scores[v.residue] = max(scores[v.residue], v.score * w_am + scores[v.residue] * (1 - w_am))
        end
    end

    return clamp.(scores, 0.0, 1.0)
end
