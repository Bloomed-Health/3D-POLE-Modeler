"""
Bayesian ACMG classification using likelihood ratio framework (Mur et al. 2023).

Implements the quantitative Bayesian framework for variant classification
where the posterior probability of pathogenicity is computed by multiplying
prior odds by the product of likelihood ratios for each evidence category.
"""

# Prior probability of pathogenicity for missense variants in POLE exonuclease domain
const PRIOR_PATHOGENICITY = 0.10  # ~10% of missense variants are pathogenic

# Likelihood ratios for ACMG evidence categories (Tavtigian et al. 2018 / Mur 2023)
const EVIDENCE_LR = Dict{Symbol, Float64}(
    # Pathogenic evidence (LR > 1)
    :PVS1 => 350.0,    # Very strong pathogenic (null variant)
    :PS1  => 18.7,     # Strong pathogenic (same AA change known pathogenic)
    :PS2  => 18.7,     # Strong pathogenic (de novo confirmed)
    :PS3  => 18.7,     # Strong pathogenic (functional studies)
    :PS4  => 18.7,     # Strong pathogenic (prevalence in affected)
    :PM1  => 4.3,      # Moderate pathogenic (mutational hotspot)
    :PM2  => 4.3,      # Moderate pathogenic (absent from population)
    :PM4  => 4.3,      # Moderate pathogenic (protein length change)
    :PM5  => 4.3,      # Moderate pathogenic (novel missense at known position)
    :PP1  => 2.08,     # Supporting pathogenic (co-segregation)
    :PP2  => 2.08,     # Supporting pathogenic (missense in low-rate gene)
    :PP3  => 2.08,     # Supporting pathogenic (computational evidence)
    :PP4  => 2.08,     # Supporting pathogenic (phenotype specific)

    # Benign evidence (LR < 1, stored as 1/LR for convenience)
    :BA1  => 1/350.0,  # Stand-alone benign (MAF > 5%)
    :BS1  => 1/18.7,   # Strong benign (MAF too high)
    :BS2  => 1/18.7,   # Strong benign (observed in healthy)
    :BS3  => 1/18.7,   # Strong benign (functional studies)
    :BS4  => 1/18.7,   # Strong benign (lack of segregation)
    :BP1  => 1/2.08,   # Supporting benign (missense in truncation gene)
    :BP3  => 1/2.08,   # Supporting benign (in-frame in non-functional)
    :BP4  => 1/2.08,   # Supporting benign (computational benign)
    :BP7  => 1/2.08,   # Supporting benign (synonymous, no splice)
)

# ACMG classification thresholds (posterior probability)
const ACMG_THRESHOLDS = (
    pathogenic = 0.99,
    likely_pathogenic = 0.90,
    likely_benign = 0.10,
    benign = 0.01,
)

"""
    ACMGClassification

Result of Bayesian variant classification.
"""
struct ACMGClassification
    acmg_class::String           # "Pathogenic", "Likely_pathogenic", "VUS", "Likely_benign", "Benign"
    posterior::Float64           # P(pathogenic | evidence)
    ci_lower::Float64           # 95% CI lower bound
    ci_upper::Float64           # 95% CI upper bound
    prior::Float64              # Prior P(pathogenic)
    evidence::Vector{Symbol}    # Evidence codes applied
    lr_product::Float64         # Product of all likelihood ratios
end

"""
    classify_variant(evidence; prior=PRIOR_PATHOGENICITY) -> ACMGClassification

Classify a variant using Bayesian framework.
Multiplies prior odds by product of likelihood ratios to get posterior.

# Arguments
- `evidence::Vector{Symbol}`: ACMG evidence codes (e.g., [:PS3, :PM1, :PP3])
- `prior::Float64`: Prior probability of pathogenicity
"""
function classify_variant(evidence::Vector{Symbol}; prior::Float64=PRIOR_PATHOGENICITY)
    # Prior odds
    prior_odds = prior / (1.0 - prior)

    # Product of likelihood ratios
    lr_product = 1.0
    for e in evidence
        lr = get(EVIDENCE_LR, e, 1.0)
        lr_product *= lr
    end

    # Posterior odds
    posterior_odds = prior_odds * lr_product

    # Convert back to probability
    posterior = posterior_odds / (1.0 + posterior_odds)
    posterior = clamp(posterior, 0.0, 1.0)

    # Approximate 95% CI using log-normal approximation
    # Variance of log(LR) product scales with number of evidence items
    n_evidence = length(evidence)
    log_var = n_evidence * 0.5  # approximate variance per evidence item
    log_posterior_odds = log(posterior_odds)
    ci_lower_odds = exp(log_posterior_odds - 1.96 * sqrt(log_var))
    ci_upper_odds = exp(log_posterior_odds + 1.96 * sqrt(log_var))
    ci_lower = ci_lower_odds / (1.0 + ci_lower_odds)
    ci_upper = ci_upper_odds / (1.0 + ci_upper_odds)
    ci_lower = clamp(ci_lower, 0.0, 1.0)
    ci_upper = clamp(ci_upper, 0.0, 1.0)

    # Determine ACMG class from posterior
    acmg_class = posterior_to_class(posterior)

    return ACMGClassification(
        acmg_class, round(posterior, digits=4),
        round(ci_lower, digits=4), round(ci_upper, digits=4),
        prior, evidence, round(lr_product, digits=2),
    )
end

"""
    posterior_to_class(posterior) -> String

Map posterior probability to ACMG classification.
"""
function posterior_to_class(posterior::Float64)
    if posterior >= ACMG_THRESHOLDS.pathogenic
        return "Pathogenic"
    elseif posterior >= ACMG_THRESHOLDS.likely_pathogenic
        return "Likely_pathogenic"
    elseif posterior <= ACMG_THRESHOLDS.benign
        return "Benign"
    elseif posterior <= ACMG_THRESHOLDS.likely_benign
        return "Likely_benign"
    else
        return "VUS"
    end
end

"""
    get_pole_evidence(mutation_id) -> Vector{Symbol}

Get pre-curated evidence codes for known POLE mutations.
Based on published functional studies and clinical data.
"""
function get_pole_evidence(mutation_id::AbstractString)
    evidence_map = Dict(
        "P286R" => [:PS3, :PS4, :PM1, :PM2, :PP3],  # Functional + clinical + hotspot
        "V411L" => [:PS3, :PS4, :PM1, :PM2, :PP3],  # Functional + clinical + hotspot
        "S297F" => [:PS3, :PM1, :PM2, :PP3, :PP4],  # Functional + PPAP phenotype
        "D275A" => [:PS3, :PM1, :PM2, :PP3],         # Functional (engineered)
        "D368A" => [:PS3, :PM1, :PM2, :PP3],         # Functional (engineered)
        "F367S" => [:PM1, :PM2, :PM5, :PP3, :PP4],  # Adjacent to catalytic, PPAP
        "L424V" => [:PS3, :PM1, :PM2, :PP3],         # Moderate functional effect
        "M444K" => [:PM1, :PM2, :PP3, :PP4],         # PPAP, computational
    )
    return get(evidence_map, mutation_id, Symbol[])
end

"""
    classify_pole_mutations(mutations) -> Vector{ACMGClassification}

Classify all known POLE mutations using pre-curated evidence.
"""
function classify_pole_mutations(mutations)
    classifications = ACMGClassification[]
    for mut in mutations
        evidence = get_pole_evidence(mut.id)
        classif = classify_variant(evidence)
        push!(classifications, classif)
    end
    return classifications
end
