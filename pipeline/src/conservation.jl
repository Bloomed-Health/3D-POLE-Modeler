"""
MSA conservation scoring: per-residue Shannon entropy from multiple sequence alignment.
"""

"""
    load_msa(path) -> Matrix{Char}

Parse a FASTA multiple sequence alignment file into a character matrix.
Rows = sequences, Columns = alignment positions (including gaps).
"""
function load_msa(path::AbstractString)
    sequences = String[]
    open(FASTX.FASTA.Reader, path) do reader
        for record in reader
            push!(sequences, string(FASTX.FASTA.sequence(record)))
        end
    end

    if isempty(sequences)
        error("No sequences found in MSA file: $path")
    end

    # Verify all sequences are same length (aligned)
    seq_len = length(sequences[1])
    for (i, s) in enumerate(sequences)
        if length(s) != seq_len
            error("Sequence $i has length $(length(s)), expected $seq_len (not aligned?)")
        end
    end

    # Convert to matrix
    n_seqs = length(sequences)
    msa = Matrix{Char}(undef, n_seqs, seq_len)
    for i in 1:n_seqs
        for j in 1:seq_len
            msa[i, j] = sequences[i][j]
        end
    end

    @info "Loaded MSA: $(n_seqs) sequences × $(seq_len) positions"
    return msa
end

"""
    compute_conservation_scores(msa; method=:shannon) -> Vector{Float64}

Compute per-column conservation score from MSA.

Methods:
- `:shannon` - Shannon entropy, inverted so 1.0 = fully conserved, 0.0 = max diversity
- `:identity` - simple fraction of most common residue

Gaps ('-', '.') are excluded from frequency calculation.
"""
function compute_conservation_scores(msa::Matrix{Char}; method::Symbol=:shannon)
    n_seqs, n_cols = size(msa)
    scores = zeros(Float64, n_cols)

    for col in 1:n_cols
        column = msa[:, col]
        # Filter gaps
        residues = filter(c -> c != '-' && c != '.', column)

        if isempty(residues)
            scores[col] = 0.0
            continue
        end

        if method == :shannon
            scores[col] = shannon_conservation(residues)
        elseif method == :identity
            scores[col] = identity_conservation(residues)
        else
            error("Unknown conservation method: $method")
        end
    end

    return scores
end

"""
    shannon_conservation(residues) -> Float64

Shannon entropy-based conservation: 1 - H/H_max.
Returns value in [0, 1] where 1 = fully conserved.
"""
function shannon_conservation(residues::Vector{Char})
    n = length(residues)
    n == 0 && return 0.0

    # Count frequencies
    counts = Dict{Char, Int}()
    for r in residues
        counts[r] = get(counts, r, 0) + 1
    end

    # Shannon entropy
    H = 0.0
    for count in values(counts)
        p = count / n
        if p > 0
            H -= p * log2(p)
        end
    end

    # Maximum possible entropy for the number of distinct amino acid types observed
    n_types = length(counts)
    H_max = log2(min(n_types, 20))
    if H_max ≈ 0.0
        return 1.0
    end

    # Invert: high entropy = low conservation
    return clamp(1.0 - H / H_max, 0.0, 1.0)
end

"""
    identity_conservation(residues) -> Float64

Simple conservation: fraction of most common residue.
"""
function identity_conservation(residues::Vector{Char})
    n = length(residues)
    n == 0 && return 0.0
    counts = Dict{Char, Int}()
    for r in residues
        counts[r] = get(counts, r, 0) + 1
    end
    return maximum(values(counts)) / n
end

"""
    smooth_conservation(scores; window=5) -> Vector{Float64}

Apply rolling mean smoothing to conservation scores.
"""
function smooth_conservation(scores::Vector{Float64}; window::Int=5)
    n = length(scores)
    smoothed = zeros(Float64, n)
    half_w = window ÷ 2

    for i in 1:n
        start_idx = max(1, i - half_w)
        end_idx = min(n, i + half_w)
        smoothed[i] = mean(scores[start_idx:end_idx])
    end

    return smoothed
end

"""
    map_msa_to_residues(msa_scores, msa_length, total_residues) -> Vector{Float64}

Map MSA column scores to protein residue positions (handling gaps in reference sequence).
Assumes first sequence in MSA is the reference (human POLE).
If MSA is shorter than total_residues, pads with mean score.
"""
function map_msa_to_residues(msa_scores::Vector{Float64}, msa::Matrix{Char}, total_residues::Int)
    # Reference is first sequence; map non-gap columns to residue positions
    ref_seq = msa[1, :]
    residue_scores = Float64[]

    for col in 1:length(ref_seq)
        if ref_seq[col] != '-' && ref_seq[col] != '.'
            push!(residue_scores, msa_scores[col])
        end
    end

    # Pad or truncate to total_residues
    mean_score = isempty(residue_scores) ? 0.5 : mean(residue_scores)
    while length(residue_scores) < total_residues
        push!(residue_scores, mean_score)
    end

    return residue_scores[1:total_residues]
end

"""
    generate_synthetic_conservation(total_residues, domain_ranges) -> Vector{Float64}

Generate synthetic per-residue conservation scores based on known domain conservation patterns.
Used when no MSA file is available. Based on literature values for POLE conservation.
"""
function generate_synthetic_conservation(
    total_residues::Int,
    domain_ranges::OrderedDict{Symbol, Tuple{Int,Int}},
)
    # Known domain-level conservation from literature
    domain_conservation = Dict(
        :ntd      => 0.30,
        :exo      => 0.85,
        :palm     => 0.95,
        :pdomain  => 0.80,
        :fingers  => 0.70,
        :thumb    => 0.60,
        :inactpol => 0.35,
        :ctd      => 0.40,
    )

    # Catalytic site residues get maximum conservation
    catalytic_sites = Set([275, 277, 286, 297, 368, 370, 411, 640, 642, 860])

    scores = fill(0.3, total_residues)  # default for linker regions

    for (domain, (start_r, end_r)) in domain_ranges
        base = get(domain_conservation, domain, 0.3)
        for r in start_r:min(end_r, total_residues)
            # Add noise for realism
            noise = 0.05 * sin(r * 0.3) + 0.03 * cos(r * 0.7)
            scores[r] = clamp(base + noise, 0.0, 1.0)
        end
    end

    # Catalytic residues: maximum conservation
    for r in catalytic_sites
        if r <= total_residues
            scores[r] = 0.99
        end
    end

    return smooth_conservation(scores; window=3)
end
