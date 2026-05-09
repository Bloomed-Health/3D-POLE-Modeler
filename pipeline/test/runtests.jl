using Test
using POLEPipeline
using OrderedCollections: OrderedDict

@testset "POLEPipeline" begin

    @testset "Config" begin
        config = PipelineConfig()
        @test config.total_residues == 2286
        @test config.primary_pdb == "4M8O"
        @test length(config.domain_ranges) == 8
        @test config.domain_ranges[:exo] == (268, 471)
        @test config.domain_ranges[:palm] == (600, 870)
    end

    @testset "Domain Ranges" begin
        # All domains should be non-overlapping major regions
        @test POLEPipeline.DOMAIN_RANGES[:ntd][1] == 1
        @test POLEPipeline.DOMAIN_RANGES[:ctd][2] == 2286
        @test length(POLEPipeline.POLE_MUTATIONS) == 8
    end

    @testset "Conservation" begin
        # Shannon conservation
        # Fully conserved column
        @test POLEPipeline.shannon_conservation(['A', 'A', 'A', 'A']) ≈ 1.0
        # Maximally diverse
        diverse = collect("ACDEFGHIKLMNPQRSTVWY")
        score = POLEPipeline.shannon_conservation(diverse)
        @test score < 0.1  # Near zero conservation

        # Smoothing
        scores = [0.0, 0.0, 1.0, 0.0, 0.0]
        smoothed = POLEPipeline.smooth_conservation(scores; window=3)
        @test smoothed[3] < 1.0  # Should be smoothed down
        @test smoothed[3] > smoothed[1]  # But still highest

        # Synthetic conservation
        synthetic = POLEPipeline.generate_synthetic_conservation(
            2286, POLEPipeline.DOMAIN_RANGES)
        @test length(synthetic) == 2286
        @test all(0.0 .<= synthetic .<= 1.0)
        # Catalytic sites should be highly conserved
        @test synthetic[286] > 0.8
    end

    @testset "Electrostatics" begin
        # Known charged residues
        @test POLEPipeline.amino_acid_charge('R') ≈ 1.0   # Arginine positive
        @test POLEPipeline.amino_acid_charge('K') ≈ 1.0   # Lysine positive
        @test POLEPipeline.amino_acid_charge('D') ≈ -1.0  # Aspartate negative
        @test POLEPipeline.amino_acid_charge('E') ≈ -1.0  # Glutamate negative
        @test POLEPipeline.amino_acid_charge('A') ≈ 0.0   # Alanine neutral

        # Three-letter codes
        @test POLEPipeline.amino_acid_charge("ARG") ≈ 1.0
        @test POLEPipeline.amino_acid_charge("ASP") ≈ -1.0

        # Synthetic charge map
        charge_map = POLEPipeline.generate_synthetic_charge_map(
            2286, POLEPipeline.DOMAIN_RANGES)
        @test length(charge_map) == 2286
        @test all(-1.0 .<= charge_map .<= 1.0)
    end

    @testset "Variants" begin
        # Built-in variants
        clinvar = POLEPipeline.get_known_clinvar_variants()
        @test length(clinvar) >= 8
        @test clinvar[1].residue == 286
        @test clinvar[1].score > 0.9

        cosmic = POLEPipeline.get_known_cosmic_variants()
        @test length(cosmic) >= 5

        # Classification score mapping
        @test POLEPipeline.classification_to_score("Pathogenic") ≈ 0.95
        @test POLEPipeline.classification_to_score("Likely_pathogenic") ≈ 0.80
        @test POLEPipeline.classification_to_score("Uncertain significance") ≈ 0.50
        @test POLEPipeline.classification_to_score("Benign") ≈ 0.05

        # Variant density
        variants = [POLEPipeline.Variant(286, 'P', 'R', :clinvar, "P", 0.99, 0.0, "")]
        density = POLEPipeline.compute_variant_density(variants, 2286; window=10)
        @test length(density) == 2286
        @test density[286] > 0.0
        @test density[1000] ≈ 0.0

        # Pathogenicity map
        path_map = POLEPipeline.compute_pathogenicity_map(clinvar, cosmic, [], 2286)
        @test length(path_map) == 2286
        @test path_map[286] > 0.5  # P286R is pathogenic
    end

    @testset "Signatures" begin
        refs = POLEPipeline.load_reference_signatures()
        @test length(refs) == 4
        @test refs[1].name == :SBS10a
        @test length(refs[1].channels) == 96
        # Profiles should sum to ~1.0
        for r in refs
            @test abs(sum(r.channels) - 1.0) < 0.01
        end

        # Pre-computed attributions
        sigs = POLEPipeline.get_pole_mutation_signatures()
        @test length(sigs) == 8
        # P286R should be dominated by SBS10a/b
        @test sigs[1][:SBS10a] + sigs[1][:SBS10b] > 0.7
        # S297F should be dominated by SBS28
        @test sigs[3][:SBS28] > 0.5
    end

    @testset "Classification" begin
        # P286R with strong evidence
        evidence = [:PS3, :PS4, :PM1, :PM2, :PP3]
        classif = POLEPipeline.classify_variant(evidence)
        @test classif.acmg_class == "Pathogenic"
        @test classif.posterior > 0.99

        # No evidence = prior only
        no_evidence = Symbol[]
        classif_none = POLEPipeline.classify_variant(no_evidence)
        @test classif_none.posterior ≈ POLEPipeline.PRIOR_PATHOGENICITY atol=0.01
        @test classif_none.acmg_class == "VUS"

        # Benign evidence
        benign_evidence = [:BS3, :BP4]
        classif_benign = POLEPipeline.classify_variant(benign_evidence)
        @test classif_benign.posterior < 0.1
        @test classif_benign.acmg_class in ["Likely_benign", "Benign"]

        # Classify all POLE mutations
        classifs = POLEPipeline.classify_pole_mutations(POLEPipeline.POLE_MUTATIONS)
        @test length(classifs) == 8
        # P286R should be pathogenic
        @test classifs[1].acmg_class in ["Pathogenic", "Likely_pathogenic"]
    end

    @testset "DDG" begin
        # P286R: Pro → Arg in exo domain
        # Pro is rigid (helix-breaker), Arg is large charged → should be destabilizing
        ddg = POLEPipeline.estimate_ddg('P', 'R', :L, 0.5, 0.85)
        @test ddg > 0.0  # Destabilizing

        # Conservative substitution V → I should be mild
        ddg_mild = POLEPipeline.estimate_ddg('V', 'I', :E, 0.7, 0.3)
        @test abs(ddg_mild) < ddg  # Less destabilizing than P→R

        # D → A at catalytic site (high conservation) should be very destabilizing
        ddg_cat = POLEPipeline.estimate_ddg('D', 'A', :L, 0.6, 0.99)
        @test ddg_cat > 2.0  # Significantly destabilizing
    end

    @testset "Structure (without PDB)" begin
        # Test SS smoothing
        ss = [:H, :H, :L, :H, :H, :H, :H, :E, :E, :L]
        smoothed = POLEPipeline.smooth_ss_assignments(ss; min_length=3)
        @test smoothed[3] == :L  # Isolated L stays (it's in a non-SS run of 1)

        # B-factor normalization
        residues = [
            POLEPipeline.Residue3D(1, "ALA", 0.0, 0.0, 0.0, 10.0, :H, "A"),
            POLEPipeline.Residue3D(2, "GLY", 1.0, 0.0, 0.0, 50.0, :L, "A"),
            POLEPipeline.Residue3D(3, "VAL", 2.0, 0.0, 0.0, 30.0, :E, "A"),
        ]
        bf = POLEPipeline.compute_bfactors(residues)
        @test length(bf) == 3
        @test bf[1] ≈ 0.0   # minimum
        @test bf[2] ≈ 1.0   # maximum
        @test 0.0 < bf[3] < 1.0

        # Domain centroids
        centroids = POLEPipeline.compute_domain_centroids(residues,
            OrderedDict(:test => (1, 3)))
        @test centroids[:test] == (1.0, 0.0, 0.0)
    end

    @testset "Output JSON structure" begin
        # Test that output functions work with minimal data
        tmpdir = mktempdir()

        # Structure JSON
        POLEPipeline.write_structure_json(joinpath(tmpdir, "test_structure.json");
            domains=POLEPipeline.DOMAIN_RANGES,
            domain_names=POLEPipeline.DOMAIN_NAMES,
            domain_abbrevs=POLEPipeline.DOMAIN_ABBREVS,
            centroids=OrderedDict(:ntd => (0.0, 0.0, 0.0)),
            domain_geometry=OrderedDict{Symbol, Dict{String,Any}}(),
            domain_ss=OrderedDict{Symbol, Vector{Dict{String,Any}}}(),
            metals=POLEPipeline.MetalSite[],
            residues=POLEPipeline.Residue3D[],
            bfactor_map=Dict(:H => 0.2, :E => 0.5, :L => 0.9),
        )
        @test isfile(joinpath(tmpdir, "test_structure.json"))

        # Scores JSON
        n = 100  # small test
        POLEPipeline.write_scores_json(joinpath(tmpdir, "test_scores.json");
            bfactors=fill(0.5, n),
            conservation=fill(0.5, n),
            variant_density=fill(0.1, n),
            pathogenicity=fill(0.1, n),
            charge=fill(0.0, n),
        )
        @test isfile(joinpath(tmpdir, "test_scores.json"))

        # Mutations JSON
        sigs = POLEPipeline.get_pole_mutation_signatures()
        classifs = POLEPipeline.classify_pole_mutations(POLEPipeline.POLE_MUTATIONS)
        ddgs = fill(2.0, length(POLEPipeline.POLE_MUTATIONS))
        POLEPipeline.write_mutations_json(joinpath(tmpdir, "test_mutations.json");
            mutations=POLEPipeline.POLE_MUTATIONS,
            signatures=sigs,
            classifications=classifs,
            ddg_values=ddgs,
        )
        @test isfile(joinpath(tmpdir, "test_mutations.json"))

        # Meta JSON
        config = PipelineConfig(; output_dir=tmpdir)
        POLEPipeline.write_meta_json(joinpath(tmpdir, "test_meta.json"), config)
        @test isfile(joinpath(tmpdir, "test_meta.json"))

        # Cleanup
        rm(tmpdir; recursive=true)
    end
end
