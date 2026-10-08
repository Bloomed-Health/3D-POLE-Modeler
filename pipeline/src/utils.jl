"""Shared utility functions for the POLE pipeline."""

using Statistics: mean

"""Sliding-window smoothing with configurable aggregation function."""
function sliding_window(values::AbstractVector{<:Real}, window::Int; op=mean)
    n = length(values)
    result = similar(values, Float64)
    half_w = window ÷ 2
    for i in 1:n
        start_idx = max(1, i - half_w)
        end_idx = min(n, i + half_w)
        result[i] = op(@view values[start_idx:end_idx])
    end
    return result
end

"""Filter residues belonging to a domain range."""
function domain_residues(residues, start_r::Int, end_r::Int)
    filtered = filter(r -> start_r <= r.resnum <= end_r, residues)
    sort!(filtered, by=r -> r.resnum)
    return filtered
end

"""Squared Euclidean distance between two residue positions."""
function dist_sq(a, b)
    return (a.x - b.x)^2 + (a.y - b.y)^2 + (a.z - b.z)^2
end
