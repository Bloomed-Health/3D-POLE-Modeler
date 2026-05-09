"""
JSON serialization: write structure, scores, mutations, and meta JSON files.
"""

"""
    write_structure_json(path, domains, metals, domain_geometry, domain_ss, centroids, residues)

Write pole_structure.json with domain geometry, SS elements, metals, backbone coords.
"""
function write_structure_json(
    path::AbstractString;
    domains::OrderedDict{Symbol, Tuple{Int,Int}},
    domain_names::OrderedDict{Symbol, String},
    domain_abbrevs::OrderedDict{Symbol, String},
    centroids::OrderedDict{Symbol, NTuple{3,Float64}},
    domain_geometry::OrderedDict{Symbol, Dict{String,Any}},
    domain_ss::OrderedDict{Symbol, Vector{Dict{String,Any}}},
    metals::Vector{MetalSite},
    residues::Vector{Residue3D},
    bfactor_map::Dict{Symbol, Float64},
)
    mkpath(dirname(path))

    # Build domain data
    domain_data = OrderedDict{String, Any}()
    for (key, (start_r, end_r)) in domains
        k = string(key)
        center = collect(get(centroids, key, (0.0, 0.0, 0.0)))
        geom = get(domain_geometry, key, Dict{String,Any}())
        ss = get(domain_ss, key, Dict{String,Any}[])

        domain_data[k] = OrderedDict(
            "name" => get(domain_names, key, k),
            "abbrev" => get(domain_abbrevs, key, uppercase(k)),
            "range" => [start_r, end_r],
            "center" => center,
            "geometry" => geom,
            "ss" => ss,
        )
    end

    # Build metal sites
    metal_data = [OrderedDict(
        "element" => m.element,
        "position" => round.([m.x, m.y, m.z], digits=2),
        "coordinating_residues" => m.coordinating_residues,
        "coordination_distances" => m.coordination_distances,
    ) for m in metals]

    # Build backbone (C-alpha coords, sampled for reasonable file size)
    backbone_data = [OrderedDict(
        "resnum" => r.resnum,
        "x" => round(r.x, digits=2),
        "y" => round(r.y, digits=2),
        "z" => round(r.z, digits=2),
        "ss" => string(r.ss_type),
    ) for r in residues]

    # B-factor map for backward compat
    bfactor_compat = OrderedDict(
        "H" => get(bfactor_map, :H, 0.2),
        "E" => get(bfactor_map, :E, 0.5),
        "L" => get(bfactor_map, :L, 0.9),
    )

    output = OrderedDict(
        "version" => "1.0.0",
        "protein" => OrderedDict(
            "uniprot" => "Q07864",
            "gene" => "POLE",
            "total_residues" => TOTAL_RESIDUES,
        ),
        "domains" => domain_data,
        "domain_order" => [string(k) for k in keys(domains)],
        "metals" => metal_data,
        "backbone" => backbone_data,
        "bfactor_map" => bfactor_compat,
    )

    open(path, "w") do io
        JSON3.pretty(io, output; allow_inf=false)
    end

    fsize = filesize(path)
    @info "Wrote $path ($(round(fsize/1024, digits=1)) KB)"
    return path
end

"""
    write_scores_json(path; bfactors, conservation, variant_density, pathogenicity, charge, domain_summaries)

Write pole_scores.json with per-residue scoring arrays.
"""
function write_scores_json(
    path::AbstractString;
    bfactors::Vector{Float64},
    conservation::Vector{Float64},
    variant_density::Vector{Float64},
    pathogenicity::Vector{Float64},
    charge::Vector{Float64},
    domain_conservation::Dict{Symbol, Float64}=Dict{Symbol,Float64}(),
    domain_density::Dict{Symbol, Float64}=Dict{Symbol,Float64}(),
    domain_charge::Dict{Symbol, Float64}=Dict{Symbol,Float64}(),
)
    mkpath(dirname(path))

    # Round arrays for file size
    round3(x) = round(x, digits=3)

    # Domain-level summaries for backward compat
    compat_conservation = OrderedDict(string(k) => round(v, digits=2) for (k, v) in domain_conservation)
    compat_density = OrderedDict(string(k) => round(v, digits=2) for (k, v) in domain_density)
    compat_charge = OrderedDict(string(k) => round(v, digits=2) for (k, v) in domain_charge)

    output = OrderedDict(
        "version" => "1.0.0",
        "n_residues" => length(bfactors),
        "per_residue" => OrderedDict(
            "bfactor" => round3.(bfactors),
            "conservation" => round3.(conservation),
            "variant_density" => round3.(variant_density),
            "pathogenicity" => round3.(pathogenicity),
            "electrostatic_charge" => round3.(charge),
        ),
        "per_domain" => OrderedDict(
            "conservation" => compat_conservation,
            "variant_density" => compat_density,
            "electrostatic_charge" => compat_charge,
        ),
    )

    open(path, "w") do io
        JSON3.pretty(io, output; allow_inf=false)
    end

    fsize = filesize(path)
    @info "Wrote $path ($(round(fsize/1024, digits=1)) KB)"
    return path
end

"""
    write_mutations_json(path; mutations, signatures, classifications, ddg_values, frameshift)

Write pole_mutations.json with mutation annotations, signature attributions, ACMG classifications.
"""
function write_mutations_json(
    path::AbstractString;
    mutations,
    signatures::Vector{Dict{Symbol, Float64}},
    classifications::Vector{ACMGClassification},
    ddg_values::Vector{Float64},
    frameshift=FRAMESHIFT_VARIANT,
)
    mkpath(dirname(path))

    mutation_data = []
    for (i, mut) in enumerate(mutations)
        sig = i <= length(signatures) ? signatures[i] : Dict{Symbol,Float64}()
        classif = i <= length(classifications) ? classifications[i] : nothing
        ddg = i <= length(ddg_values) ? ddg_values[i] : 0.0

        mut_entry = OrderedDict(
            "id" => mut.id,
            "residue" => mut.residue,
            "domain" => string(mut.domain),
            "ref_aa" => string(mut.ref_aa),
            "alt_aa" => string(mut.alt_aa),
            "label" => mut.label,
            "detail" => mut.detail,
            "ddg_kcal_mol" => ddg,
            "signature_attribution" => OrderedDict(string(k) => v for (k, v) in sig),
        )

        if classif !== nothing
            # Get auditable citations for each evidence code
            cit_map = get_pole_evidence_citations(mut.id)
            evidence_detail = []
            for e in classif.evidence
                entry = OrderedDict("code" => string(e), "lr" => get(EVIDENCE_LR, e, 1.0))
                citation = get(cit_map, e, "")
                if !isempty(citation)
                    entry["citation"] = citation
                end
                push!(evidence_detail, entry)
            end

            mut_entry["classification"] = OrderedDict(
                "acmg_class" => classif.acmg_class,
                "posterior" => classif.posterior,
                "ci_lower" => classif.ci_lower,
                "ci_upper" => classif.ci_upper,
                "evidence" => evidence_detail,
                "lr_product" => classif.lr_product,
                "framework" => "Mur et al. 2023 (Genome Med 15:85) / Tavtigian et al. 2018",
            )
        end

        push!(mutation_data, mut_entry)
    end

    # Frameshift data
    fs_data = OrderedDict(
        "hgvs_c" => frameshift.hgvs_c,
        "hgvs_p" => frameshift.hgvs_p,
        "truncation_site" => frameshift.truncation_site,
        "total_residues" => TOTAL_RESIDUES,
        "wt_codons" => [OrderedDict("codon" => c.codon, "aa" => c.aa, "pos" => c.pos)
                        for c in frameshift.wt_codons],
        "mut_codons" => [OrderedDict("codon" => c.codon, "aa" => c.aa, "pos" => c.pos, "status" => c.status)
                         for c in frameshift.mut_codons],
    )

    output = OrderedDict(
        "version" => "1.0.0",
        "mutations" => mutation_data,
        "frameshift" => fs_data,
    )

    open(path, "w") do io
        JSON3.pretty(io, output; allow_inf=false)
    end

    fsize = filesize(path)
    @info "Wrote $path ($(round(fsize/1024, digits=1)) KB)"
    return path
end

"""
    write_meta_json(path, config; pdb_resolutions, checksums)

Write pole_meta.json with pipeline provenance metadata.
"""
function write_meta_json(
    path::AbstractString,
    config::PipelineConfig;
    pdb_resolutions::Dict{String, Float64}=Dict{String,Float64}(),
    checksums::Dict{String, String}=Dict{String,String}(),
)
    mkpath(dirname(path))

    output = OrderedDict(
        "version" => "1.0.0",
        "pipeline" => OrderedDict(
            "name" => "POLEPipeline",
            "version" => "0.1.0",
            "julia_version" => string(VERSION),
            "run_timestamp" => Dates.format(now(), "yyyy-mm-ddTHH:MM:SS"),
        ),
        "data_sources" => OrderedDict(
            "primary_pdb" => config.primary_pdb,
            "all_pdb_ids" => config.pdb_ids,
            "pdb_resolutions" => pdb_resolutions,
        ),
        "parameters" => OrderedDict(
            "conservation_window" => config.conservation_window,
            "density_window" => config.density_window,
            "charge_window" => config.charge_window,
            "prior_pathogenicity" => PRIOR_PATHOGENICITY,
        ),
        "checksums" => checksums,
    )

    open(path, "w") do io
        JSON3.pretty(io, output; allow_inf=false)
    end

    @info "Wrote $path"
    return path
end

"""
    compute_file_checksum(path) -> String

Compute SHA-256 checksum of a file.
"""
function compute_file_checksum(path::AbstractString)
    return bytes2hex(open(SHA.sha256, path))
end

"""
    compute_domain_summary(scores, domain_ranges) -> Dict{Symbol, Float64}

Compute per-domain mean score from per-residue array.
"""
function compute_domain_summary(
    scores::Vector{Float64},
    domain_ranges::OrderedDict{Symbol, Tuple{Int,Int}},
)
    result = Dict{Symbol, Float64}()
    for (domain, (start_r, end_r)) in domain_ranges
        end_idx = min(end_r, length(scores))
        start_idx = max(1, start_r)
        if start_idx <= end_idx
            result[domain] = mean(scores[start_idx:end_idx])
        else
            result[domain] = 0.0
        end
    end
    return result
end
