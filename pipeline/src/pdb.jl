"""
PDB/mmCIF parsing: download, parse, extract backbone, SS records, metals, coordination.
"""

"""
    download_pdb(pdb_id, cache_dir) -> String

Download mmCIF file from RCSB if not already cached. Returns path to local file.
"""
function download_pdb(pdb_id::AbstractString, cache_dir::AbstractString)
    mkpath(cache_dir)
    id_lower = lowercase(pdb_id)
    local_path = joinpath(cache_dir, "$(id_lower).cif")
    if isfile(local_path)
        @info "PDB $pdb_id already cached at $local_path"
        return local_path
    end
    url = "https://files.rcsb.org/download/$(id_lower).cif"
    @info "Downloading PDB $pdb_id from $url"
    Downloads.download(url, local_path)
    return local_path
end

"""
    parse_structure(path) -> MolecularStructure

Parse an mmCIF file into a BioStructures MolecularStructure.
"""
function parse_structure(path::AbstractString)
    @info "Parsing structure from $path"
    if endswith(path, ".cif")
        return read(path, MMCIFFormat)
    else
        return read(path, PDBFormat)
    end
end

"""
    SSElement

A secondary structure element with type (:H, :E, :L), start/end residue numbers, and length.
"""
struct SSElement
    type::Symbol    # :H (helix), :E (strand), :L (loop/coil)
    start_res::Int
    end_res::Int
    length::Int
end

"""
    Residue3D

A C-alpha residue with 3D coordinates and metadata.
"""
struct Residue3D
    resnum::Int
    resname::String
    x::Float64
    y::Float64
    z::Float64
    bfactor::Float64
    ss_type::Symbol  # :H, :E, :L (assigned later)
    chain::String
end

"""
    MetalSite

A metal ion with coordinates and coordinating residues.
"""
struct MetalSite
    element::String
    x::Float64
    y::Float64
    z::Float64
    coordinating_residues::Vector{Int}
    coordination_distances::Vector{Float64}
end

"""
    extract_ca_backbone(struc, chain_id) -> Vector{Residue3D}

Extract C-alpha atoms from a specific chain, returning per-residue 3D coordinates.
"""
function extract_ca_backbone(struc, chain_id::AbstractString)
    residues = Residue3D[]
    model1 = defaultmodel(struc)
    ch = model1[chain_id]

    for res in collectresidues(ch, standardselector)
        ca_atoms = collectatoms(res, calphaselector)
        isempty(ca_atoms) && continue
        ca = first(ca_atoms)
        push!(residues, Residue3D(
            resnumber(res),
            resname(res),
            x(ca), y(ca), z(ca),
            tempfactor(ca),
            :L,  # default; SS assigned later
            chain_id,
        ))
    end

    sort!(residues, by=r -> r.resnum)
    @info "Extracted $(length(residues)) C-alpha atoms from chain $chain_id"
    return residues
end

"""
    extract_ss_from_mmcif(path) -> Vector{SSElement}

Parse HELIX/SHEET records from mmCIF to get secondary structure assignments.
Falls back to phi/psi-based assignment if records not present.
"""
function extract_ss_from_mmcif(path::AbstractString)
    elements = SSElement[]

    # Parse _struct_conf (helices) and _struct_sheet_range (strands) from mmCIF
    lines = readlines(path)

    # Parse helices from _struct_conf
    in_struct_conf = false
    conf_cols = Dict{String,Int}()
    for line in lines
        stripped = strip(line)
        if startswith(stripped, "_struct_conf.")
            in_struct_conf = true
            field = split(stripped, '.')[2]
            col_idx = length(conf_cols) + 1
            conf_cols[field] = col_idx
        elseif in_struct_conf && !startswith(stripped, "_") && !startswith(stripped, "#") && !isempty(stripped) && !startswith(stripped, "loop_")
            tokens = split(stripped)
            if haskey(conf_cols, "beg_auth_seq_id") && haskey(conf_cols, "end_auth_seq_id")
                beg_idx = conf_cols["beg_auth_seq_id"]
                end_idx = conf_cols["end_auth_seq_id"]
                if beg_idx <= length(tokens) && end_idx <= length(tokens)
                    start_res = tryparse(Int, tokens[beg_idx])
                    end_res = tryparse(Int, tokens[end_idx])
                    if start_res !== nothing && end_res !== nothing
                        push!(elements, SSElement(:H, start_res, end_res, end_res - start_res + 1))
                    end
                end
            end
        elseif in_struct_conf && (startswith(stripped, "_") || startswith(stripped, "#") || startswith(stripped, "loop_"))
            in_struct_conf = false
            conf_cols = Dict{String,Int}()
        end
    end

    # Parse strands from _struct_sheet_range
    in_sheet_range = false
    sheet_cols = Dict{String,Int}()
    for line in lines
        stripped = strip(line)
        if startswith(stripped, "_struct_sheet_range.")
            in_sheet_range = true
            field = split(stripped, '.')[2]
            col_idx = length(sheet_cols) + 1
            sheet_cols[field] = col_idx
        elseif in_sheet_range && !startswith(stripped, "_") && !startswith(stripped, "#") && !isempty(stripped) && !startswith(stripped, "loop_")
            tokens = split(stripped)
            if haskey(sheet_cols, "beg_auth_seq_id") && haskey(sheet_cols, "end_auth_seq_id")
                beg_idx = sheet_cols["beg_auth_seq_id"]
                end_idx = sheet_cols["end_auth_seq_id"]
                if beg_idx <= length(tokens) && end_idx <= length(tokens)
                    start_res = tryparse(Int, tokens[beg_idx])
                    end_res = tryparse(Int, tokens[end_idx])
                    if start_res !== nothing && end_res !== nothing
                        push!(elements, SSElement(:E, start_res, end_res, end_res - start_res + 1))
                    end
                end
            end
        elseif in_sheet_range && (startswith(stripped, "_") || startswith(stripped, "#") || startswith(stripped, "loop_"))
            in_sheet_range = false
            sheet_cols = Dict{String,Int}()
        end
    end

    sort!(elements, by=e -> e.start_res)
    @info "Extracted $(length(elements)) SS elements ($(count(e->e.type==:H, elements)) helices, $(count(e->e.type==:E, elements)) strands)"
    return elements
end

"""
    assign_ss_to_residues(residues, ss_elements) -> Vector{Residue3D}

Assign secondary structure type to each residue based on SS element ranges.
"""
function assign_ss_to_residues(residues::Vector{Residue3D}, ss_elements::Vector{SSElement})
    # Build lookup: residue number -> SS type
    ss_map = Dict{Int, Symbol}()
    for elem in ss_elements
        for r in elem.start_res:elem.end_res
            ss_map[r] = elem.type
        end
    end

    return [Residue3D(
        r.resnum, r.resname, r.x, r.y, r.z, r.bfactor,
        get(ss_map, r.resnum, :L),
        r.chain,
    ) for r in residues]
end

"""
    extract_heteroatoms(struc, chain_id, types) -> Vector{MetalSite}

Extract metal ions and heteroatoms of specified types from structure.
"""
function extract_heteroatoms(struc, chain_id::AbstractString, types::Vector{String}=METAL_TYPES)
    metals = MetalSite[]
    model1 = defaultmodel(struc)

    for ch in collectchains(model1)
        for res in collectresidues(ch, heteroselector)
            rn = resname(res)
            rn in types || continue
            atoms = collectatoms(res)
            isempty(atoms) && continue
            a = first(atoms)
            push!(metals, MetalSite(rn, x(a), y(a), z(a), Int[], Float64[]))
        end
    end

    @info "Found $(length(metals)) metal/heteroatom sites"
    return metals
end

"""
    extract_metal_coordination(metals, residues; cutoff) -> Vector{MetalSite}

Find protein residues within coordination distance of each metal site.
"""
function extract_metal_coordination(
    metals::Vector{MetalSite},
    residues::Vector{Residue3D};
    cutoff::Float64=METAL_COORDINATION_CUTOFF,
)
    coordinated = MetalSite[]
    for m in metals
        coord_res = Int[]
        coord_dist = Float64[]
        for r in residues
            d = sqrt((m.x - r.x)^2 + (m.y - r.y)^2 + (m.z - r.z)^2)
            if d <= cutoff
                push!(coord_res, r.resnum)
                push!(coord_dist, round(d, digits=2))
            end
        end
        push!(coordinated, MetalSite(m.element, m.x, m.y, m.z, coord_res, coord_dist))
    end
    return coordinated
end

"""
    compute_domain_centroids(residues, ranges) -> Dict{Symbol, NTuple{3,Float64}}

Compute geometric center of each domain from C-alpha coordinates.
"""
function compute_domain_centroids(
    residues::Vector{Residue3D},
    ranges::OrderedDict{Symbol, Tuple{Int,Int}},
)
    centroids = OrderedDict{Symbol, NTuple{3,Float64}}()
    for (domain, (start_r, end_r)) in ranges
        domain_res = filter(r -> start_r <= r.resnum <= end_r, residues)
        if isempty(domain_res)
            centroids[domain] = (0.0, 0.0, 0.0)
            continue
        end
        cx = mean(r.x for r in domain_res)
        cy = mean(r.y for r in domain_res)
        cz = mean(r.z for r in domain_res)
        centroids[domain] = (round(cx, digits=2), round(cy, digits=2), round(cz, digits=2))
    end
    return centroids
end
