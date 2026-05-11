#!/usr/bin/env julia
"""
POLEPipeline entry point.
Run with: julia --project=pipeline pipeline/run.jl

Produces 4 JSON files in data/:
  - pole_structure.json  (domains, SS elements, metals, backbone)
  - pole_scores.json     (per-residue bfactor, conservation, density, pathogenicity, charge)
  - pole_mutations.json  (mutations, signatures, ACMG classification, DDG)
  - pole_meta.json       (provenance, checksums, timestamps)
"""

using Pkg
Pkg.activate(joinpath(@__DIR__))

using POLEPipeline
using OrderedCollections: OrderedDict

function run_pipeline(; use_pdb::Bool=true)
    config = PipelineConfig()

    @info "=== POLEPipeline v0.1.0 ==="
    @info "Output directory: $(config.output_dir)"
    mkpath(config.output_dir)
    mkpath(config.cache_dir)

    # ═══════════════════════════════════════════════════════════════════
    # Phase 1: Structure
    # ═══════════════════════════════════════════════════════════════════
    @info "─── Phase 1: Structure ───"

    residues = POLEPipeline.Residue3D[]
    ss_elements = POLEPipeline.SSElement[]
    metals = POLEPipeline.MetalSite[]

    if use_pdb
        try
            # Download and parse primary PDB structure
            pdb_path = POLEPipeline.download_pdb(config.primary_pdb, config.cache_dir)
            struc = POLEPipeline.parse_structure(pdb_path)

            # Extract C-alpha backbone
            residues = POLEPipeline.extract_ca_backbone(struc, config.primary_chain)

            # Secondary structure assignment: DSSP → mmCIF records → geometry heuristic
            ss_elements = POLEPipeline.SSElement[]

            # Try DSSP binary first (most accurate)
            dssp_result = POLEPipeline.run_dssp(pdb_path)
            if dssp_result !== nothing && !isempty(dssp_result)
                ss_elements = dssp_result
                @info "Using DSSP binary for SS assignment"
            else
                # Fallback: mmCIF _struct_conf records
                ss_elements = POLEPipeline.extract_ss_from_mmcif(pdb_path)
                if !isempty(ss_elements)
                    @info "Using mmCIF SS records for SS assignment"
                end
            end

            # Assign SS to residues
            if !isempty(ss_elements)
                residues = POLEPipeline.assign_ss_to_residues(residues, ss_elements)
            else
                # Final fallback: geometry-based assignment
                @info "Using geometry heuristic for SS assignment"
                residues = POLEPipeline.assign_secondary_structure(residues)
            end

            # Extract metals
            metals = POLEPipeline.extract_heteroatoms(struc, config.primary_chain)
            metals = POLEPipeline.extract_metal_coordination(metals, residues)

            @info "Structure loaded: $(length(residues)) residues, $(length(metals)) metal sites"
        catch e
            @warn "PDB loading failed, continuing with synthetic data" exception=e
            use_pdb = false
        end
    end

    # Compute domain centroids
    centroids = if !isempty(residues)
        POLEPipeline.compute_domain_centroids(residues, config.domain_ranges)
    else
        # Fallback centroids from existing viewer
        OrderedDict(
            :ntd      => (-18.0, 8.0, -5.0),
            :exo      => (-10.0, 2.0, 7.0),
            :palm     => (0.0, -1.0, 0.0),
            :pdomain  => (4.0, 1.0, -3.0),
            :fingers  => (7.0, 5.0, 4.0),
            :thumb    => (3.0, -9.0, 7.0),
            :inactpol => (10.0, -4.0, 1.0),
            :ctd      => (16.0, 1.0, -5.0),
        )
    end

    # Compute domain geometry and SS elements
    domain_geometry = if !isempty(residues)
        POLEPipeline.compute_domain_geometry(residues, config.domain_ranges)
    else
        OrderedDict{Symbol, Dict{String,Any}}()
    end

    domain_ss = if !isempty(residues)
        POLEPipeline.compute_domain_ss_elements(residues, config.domain_ranges)
    else
        OrderedDict{Symbol, Vector{Dict{String,Any}}}()
    end

    # B-factor computation
    norm_bfactors = if !isempty(residues)
        POLEPipeline.compute_bfactors(residues)
    else
        Float64[]
    end

    bfactor_map = if !isempty(residues)
        POLEPipeline.compute_per_ss_bfactor(residues, norm_bfactors)
    else
        Dict(:H => 0.2, :E => 0.5, :L => 0.9)
    end

    # Write structure JSON
    structure_path = joinpath(config.output_dir, "pole_structure.json")
    POLEPipeline.write_structure_json(structure_path;
        domains=config.domain_ranges,
        domain_names=POLEPipeline.DOMAIN_NAMES,
        domain_abbrevs=POLEPipeline.DOMAIN_ABBREVS,
        centroids=centroids,
        domain_geometry=domain_geometry,
        domain_ss=domain_ss,
        metals=metals,
        residues=residues,
        bfactor_map=bfactor_map,
    )

    # ═══════════════════════════════════════════════════════════════════
    # Phase 1b: AlphaFold
    # ═══════════════════════════════════════════════════════════════════
    @info "─── Phase 1b: AlphaFold ───"

    plddt = Float64[]
    try
        af_path = POLEPipeline.download_alphafold(
            POLEPipeline.ALPHAFOLD_UNIPROT, config.cache_dir)
        af_struc = POLEPipeline.parse_structure(af_path)
        plddt = POLEPipeline.extract_plddt_scores(af_struc, "A";
            total_residues=config.total_residues)
        @info "AlphaFold pLDDT: mean=$(round(mean(plddt), digits=2)), min=$(round(minimum(plddt), digits=2)), max=$(round(maximum(plddt), digits=2))"
    catch e
        @warn "AlphaFold download/parsing failed, continuing without pLDDT" exception=e
    end

    # ═══════════════════════════════════════════════════════════════════
    # Phase 2: Per-Residue Scores
    # ═══════════════════════════════════════════════════════════════════
    @info "─── Phase 2: Per-Residue Scores ───"

    # Conservation
    msa_path = joinpath(config.data_dir, "msa", "pole_orthologs.fasta")
    conservation = if isfile(msa_path)
        msa = POLEPipeline.load_msa(msa_path)
        raw_scores = POLEPipeline.compute_conservation_scores(msa)
        mapped = POLEPipeline.map_msa_to_residues(raw_scores, msa, config.total_residues)
        POLEPipeline.smooth_conservation(mapped; window=config.conservation_window)
    else
        @info "No MSA file found, generating synthetic conservation scores"
        POLEPipeline.generate_synthetic_conservation(config.total_residues, config.domain_ranges)
    end

    # Variant loading
    clinvar_path = joinpath(config.data_dir, "variants", "clinvar_pole.tsv")
    cosmic_path = joinpath(config.data_dir, "variants", "cosmic_pole.tsv")
    am_path = joinpath(config.data_dir, "variants", "alphamissense_pole.tsv")

    clinvar = POLEPipeline.load_clinvar_variants(clinvar_path)
    cosmic = POLEPipeline.load_cosmic_variants(cosmic_path)
    alphamissense = POLEPipeline.load_alphamissense_scores(am_path)

    # Variant density
    all_variants = vcat(clinvar, cosmic, alphamissense)
    variant_density = POLEPipeline.compute_variant_density(
        all_variants, config.total_residues; window=config.density_window)

    # Pathogenicity map
    pathogenicity = POLEPipeline.compute_pathogenicity_map(
        clinvar, cosmic, alphamissense, config.total_residues)

    # Electrostatic charge
    charge = if !isempty(residues)
        POLEPipeline.compute_charge_map(residues; window=config.charge_window)
    else
        POLEPipeline.generate_synthetic_charge_map(config.total_residues, config.domain_ranges)
    end

    # Per-residue B-factors (pad to full length if structure partial)
    full_bfactors = if !isempty(norm_bfactors)
        # Map structural bfactors to full residue range
        bf = fill(0.5, config.total_residues)
        for (i, r) in enumerate(residues)
            if 1 <= r.resnum <= config.total_residues
                bf[r.resnum] = norm_bfactors[i]
            end
        end
        bf
    else
        fill(0.5, config.total_residues)
    end

    # Domain summaries for backward compat
    domain_conservation = POLEPipeline.compute_domain_summary(conservation, config.domain_ranges)
    domain_density = POLEPipeline.compute_domain_summary(variant_density, config.domain_ranges)
    domain_charge = POLEPipeline.compute_domain_summary(charge, config.domain_ranges)

    # Write scores JSON
    scores_path = joinpath(config.output_dir, "pole_scores.json")
    POLEPipeline.write_scores_json(scores_path;
        bfactors=full_bfactors,
        conservation=conservation,
        variant_density=variant_density,
        pathogenicity=pathogenicity,
        charge=charge,
        plddt=plddt,
        domain_conservation=domain_conservation,
        domain_density=domain_density,
        domain_charge=domain_charge,
    )

    # ═══════════════════════════════════════════════════════════════════
    # Phase 3: Mutation Analysis
    # ═══════════════════════════════════════════════════════════════════
    @info "─── Phase 3: Mutation Analysis ───"

    mutations = POLEPipeline.POLE_MUTATIONS

    # Signature deconvolution
    # Prefer official COSMIC v3.4 file when present; fall back to built-in approximations
    cosmic_sigs_path = joinpath(config.data_dir, "COSMIC_v3.4_SBS_GRCh38.txt")
    sig_refs = POLEPipeline.load_reference_signatures(
        isfile(cosmic_sigs_path) ? cosmic_sigs_path : "")
    signatures = POLEPipeline.get_pole_mutation_signatures()

    # ACMG classification
    classifications = POLEPipeline.classify_pole_mutations(mutations)

    # DDG estimation
    ddg_values = if !isempty(residues)
        POLEPipeline.batch_ddg(mutations, residues, conservation)
    else
        # Estimate without structure using default burial
        [POLEPipeline.estimate_ddg(
            mut.ref_aa, mut.alt_aa, :L, 0.5,
            mut.residue <= length(conservation) ? conservation[mut.residue] : 0.5,
        ) for mut in mutations]
    end

    # Write mutations JSON
    mutations_path = joinpath(config.output_dir, "pole_mutations.json")
    POLEPipeline.write_mutations_json(mutations_path;
        mutations=mutations,
        signatures=signatures,
        classifications=classifications,
        ddg_values=ddg_values,
    )

    # ═══════════════════════════════════════════════════════════════════
    # Phase 4: Meta
    # ═══════════════════════════════════════════════════════════════════
    @info "─── Phase 4: Metadata ───"

    # Compute checksums
    checksums = Dict{String, String}()
    for fname in ["pole_structure.json", "pole_scores.json", "pole_mutations.json"]
        fpath = joinpath(config.output_dir, fname)
        if isfile(fpath)
            checksums[fname] = POLEPipeline.compute_file_checksum(fpath)
        end
    end

    meta_path = joinpath(config.output_dir, "pole_meta.json")
    POLEPipeline.write_meta_json(meta_path, config; checksums=checksums)

    @info "=== Pipeline complete ==="
    @info "Output files:"
    for fname in ["pole_structure.json", "pole_scores.json", "pole_mutations.json", "pole_meta.json"]
        fpath = joinpath(config.output_dir, fname)
        if isfile(fpath)
            @info "  $fname ($(round(filesize(fpath)/1024, digits=1)) KB)"
        end
    end
end

# Run pipeline
if abspath(PROGRAM_FILE) == @__FILE__
    # Default: try PDB but fall back gracefully
    use_pdb = !("--no-pdb" in ARGS)
    run_pipeline(; use_pdb=use_pdb)
end
