"""
Structural math: SS assignment, B-factors, domain geometry.
"""

"""
    assign_secondary_structure(residues) -> Vector{Residue3D}

DSSP-inspired secondary structure assignment from C-alpha geometry.
Uses C-alpha distance patterns to identify helices and strands when PDB SS records
are not available.

- α-helix: Cα(i) to Cα(i+3) distance ~5.0-5.5 Å
- β-strand: Cα(i) to Cα(i+2) distance ~6.5-7.0 Å (extended)
"""
function assign_secondary_structure(residues::Vector{Residue3D})
    n = length(residues)
    ss_types = fill(:L, n)

    # Build residue number -> index mapping
    idx_map = Dict(r.resnum => i for (i, r) in enumerate(residues))

    for i in 1:n
        r = residues[i]

        # Check for α-helix: i to i+3 distance ≈ 5.0-5.5 Å
        if i + 3 <= n
            r3 = residues[i+3]
            d13 = sqrt((r.x - r3.x)^2 + (r.y - r3.y)^2 + (r.z - r3.z)^2)
            if 4.5 <= d13 <= 6.0
                # Also check i to i+2 for consistency
                r2 = residues[i+2]
                d12 = sqrt((r.x - r2.x)^2 + (r.y - r2.y)^2 + (r.z - r2.z)^2)
                if 4.5 <= d12 <= 6.5
                    ss_types[i] = :H
                end
            end
        end

        # Check for β-strand: extended backbone
        if i + 2 <= n && ss_types[i] != :H
            r2 = residues[i+2]
            d = sqrt((r.x - r2.x)^2 + (r.y - r2.y)^2 + (r.z - r2.z)^2)
            if 6.0 <= d <= 7.5
                ss_types[i] = :E
            end
        end
    end

    # Smooth: remove isolated SS assignments (< 3 consecutive)
    ss_types = smooth_ss_assignments(ss_types)

    return [Residue3D(
        r.resnum, r.resname, r.x, r.y, r.z, r.bfactor,
        ss_types[i], r.chain,
    ) for (i, r) in enumerate(residues)]
end

"""
    smooth_ss_assignments(ss_types; min_length=3) -> Vector{Symbol}

Remove isolated secondary structure assignments shorter than min_length.
"""
function smooth_ss_assignments(ss_types::Vector{Symbol}; min_length::Int=3)
    n = length(ss_types)
    result = copy(ss_types)
    i = 1
    while i <= n
        current = result[i]
        j = i
        while j <= n && result[j] == current
            j += 1
        end
        run_length = j - i
        if current != :L && run_length < min_length
            result[i:j-1] .= :L
        end
        i = j
    end
    return result
end

"""
    compute_bfactors(residues) -> Vector{Float64}

Normalize crystallographic B-factors to 0-1 range.
"""
function compute_bfactors(residues::Vector{Residue3D})
    bvals = [r.bfactor for r in residues]
    if isempty(bvals)
        return Float64[]
    end
    bmin, bmax = extrema(bvals)
    if bmax ≈ bmin
        return fill(0.5, length(bvals))
    end
    return [(b - bmin) / (bmax - bmin) for b in bvals]
end

"""
    compute_domain_ss_elements(residues, ranges) -> Dict{Symbol, Vector{Dict}}

Group secondary structure into per-domain element lists.
Output format matches JS: [{t:"H", n:12}, {t:"L", n:5}, ...]
"""
function compute_domain_ss_elements(
    residues::Vector{Residue3D},
    ranges::OrderedDict{Symbol, Tuple{Int,Int}},
)
    result = OrderedDict{Symbol, Vector{Dict{String,Any}}}()

    for (domain, (start_r, end_r)) in ranges
        domain_res = filter(r -> start_r <= r.resnum <= end_r, residues)
        sort!(domain_res, by=r -> r.resnum)

        elements = Dict{String,Any}[]
        isempty(domain_res) && (result[domain] = elements; continue)

        current_type = domain_res[1].ss_type
        current_count = 1

        for i in 2:length(domain_res)
            if domain_res[i].ss_type == current_type
                current_count += 1
            else
                push!(elements, Dict("t" => string(current_type), "n" => current_count))
                current_type = domain_res[i].ss_type
                current_count = 1
            end
        end
        push!(elements, Dict("t" => string(current_type), "n" => current_count))

        result[domain] = elements
    end

    return result
end

"""
    compute_domain_geometry(residues, ranges) -> Dict{Symbol, Dict}

Compute center, bounding box, and orientation for each domain.
"""
function compute_domain_geometry(
    residues::Vector{Residue3D},
    ranges::OrderedDict{Symbol, Tuple{Int,Int}},
)
    result = OrderedDict{Symbol, Dict{String,Any}}()

    for (domain, (start_r, end_r)) in ranges
        domain_res = filter(r -> start_r <= r.resnum <= end_r, residues)
        if isempty(domain_res)
            result[domain] = Dict{String,Any}(
                "center" => [0.0, 0.0, 0.0],
                "bbox_min" => [0.0, 0.0, 0.0],
                "bbox_max" => [0.0, 0.0, 0.0],
                "radius" => 0.0,
            )
            continue
        end

        xs = [r.x for r in domain_res]
        ys = [r.y for r in domain_res]
        zs = [r.z for r in domain_res]

        cx, cy, cz = mean(xs), mean(ys), mean(zs)

        # Approximate radius: mean distance from centroid
        dists = [sqrt((r.x - cx)^2 + (r.y - cy)^2 + (r.z - cz)^2) for r in domain_res]
        radius = mean(dists)

        result[domain] = Dict{String,Any}(
            "center" => round.([cx, cy, cz], digits=2),
            "bbox_min" => round.([minimum(xs), minimum(ys), minimum(zs)], digits=2),
            "bbox_max" => round.([maximum(xs), maximum(ys), maximum(zs)], digits=2),
            "radius" => round(radius, digits=2),
        )
    end

    return result
end

"""
    compute_per_residue_bfactor_map(residues) -> Dict{Symbol, Float64}

Compute average normalized B-factor per SS type (for backward compat with JS BFACTOR_MAP).
"""
function compute_per_ss_bfactor(residues::Vector{Residue3D}, norm_bfactors::Vector{Float64})
    result = Dict{Symbol, Float64}()
    for ss_type in [:H, :E, :L]
        indices = findall(r -> r.ss_type == ss_type, residues)
        if !isempty(indices)
            result[ss_type] = round(mean(norm_bfactors[indices]), digits=3)
        end
    end
    return result
end
