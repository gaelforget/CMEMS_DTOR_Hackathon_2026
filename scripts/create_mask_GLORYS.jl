### A Pluto.jl notebook ###
# v1.0.1

using Markdown
using InteractiveUtils

# ╔═╡ 41809153-52b3-4f91-8edd-79e5522bad36
begin
	import Pkg  
    Pkg.activate(joinpath(@__DIR__, ".."))
end

# ╔═╡ 447ba694-6a19-4400-b0da-79e97b8f5c51
using MeshArrays, Glob, Statistics, NCDatasets, Dates, CairoMakie, JLD2, MITgcm

# ╔═╡ 6c84e1e0-a9b1-45a1-b499-5df954d1c6ef
md""" # Create regional masks for GLORYS grid based on ECCOr6 mask field"""

# ╔═╡ f055de04-692e-40d4-8219-3b1d8a95a2cc
grpdir = get(ENV, "GROUP", "")

# ╔═╡ db7654cb-b219-4216-ae63-cd9900f662a8
begin
	datadir = joinpath(grpdir, "GLORYS/GLOBAL_MULTIYEAR_PHY_001_030/cmems_mod_glo_phy_my_0.083deg_P1M-m_202311")
	println(datadir)
end

# ╔═╡ 760143e6-e022-4285-ae1e-2a058cbcc749
md""" #### get all filenames & read some of them"""

# ╔═╡ e546615b-ea5d-4203-964b-de68e2a0d65e
begin
    files = String[]
    
    for year in 2023:2023
        year_dir = joinpath(datadir, string(year))
    
        append!(
            files,
            glob(
                "mercatorglorys12v1_gl12_mean_$(year)*.nc",
                year_dir,
            ),
        )
    end
    sort!(files)

    @show length(files)
    # println.(files)
end

# ╔═╡ 921a9e4e-f4b1-430c-9884-5619501a5337
ds = NCDataset(files[1:1], "r";
    aggdim = "time",
    deferopen = true,
)

# ╔═╡ 1a6f854d-f49e-4603-8652-f467488f5d89
md""" ## Grid - build using MeshArrays """

# ╔═╡ dedef362-5070-4901-b322-a3d2894968f7
begin
	RC = ds["depth"][:] # depth of layer center
	lon = ds["longitude"][:]
	lat = ds["latitude"][:]
	temp = ds["thetao"][:,:,:,1]
	γ, Γ_raw = glorys_grid(lon, lat)
end

# ╔═╡ 8c4e5825-f5de-4859-9b82-cbdab561d5a5
begin
	RF = [0; 0.5*(RC[1:end-1]+RC[2:end]); RC[end]+0.5*(RC[end]-RC[end-1])] 
	DRF = diff(RF)
	hFacC = read(1.0*(!ismissing).(temp),γ)
	depth = Γ_raw.XC * 0
	for t in axes(hFacC, 1), k in axes(hFacC, 2)
		depth[t] += hFacC[t, k] .* DRF[k]
	end
	Γ = merge(Γ_raw, (RC = RC, RF = RF, DRF = DRF, hFacC = hFacC, Depth=depth))
end

# ╔═╡ cd56d103-def3-438c-a37f-e28d389df823
begin
	# ECCOr6 grid
	griddir = "/nfs/micklab005/gforget/2026_datasets/ECCO4r6_output/grid"
	γ6 = GridSpec("LatLonCap", griddir, ioPrec=Float32)
	Γ6 = GridLoad(γ6;option="full")
end

# ╔═╡ f4668412-0053-44d4-8ba4-9eed43e7c47a


# ╔═╡ 3c9add6e-62ed-44f1-b8ff-ea50b9991873


# ╔═╡ e84f674a-4e24-48cd-8dbb-53c345f71ee9
md""" ## get masks for GLORYS """

# ╔═╡ 5f511f7b-4cff-4691-961b-21b619319a33
maskdir = expanduser("~/data_small")

# ╔═╡ 2f6bd8b9-1123-43a8-a769-80fb2139dbb0
begin
    # Define region specifications
    region_names = ("lab",)
    nregs = length(region_names)
    
    nx6, ny6, nz6 = size(write(Γ6.hFacC))
    
    # Read all masks efficiently
    read_mask(fname) = begin
        raw = Array{Float32}(undef, nx6, ny6)
        MITgcm.read!(fname, raw)
        raw .= ntoh.(raw)
        γ6.read(raw, Γ6.XC)
    end

    # refine_inter_maskC, etc
    region_masks_ref = [
        read_mask(joinpath(maskdir, "refine_$(name)_maskC"))
        for name in region_names
    ]
end

# ╔═╡ 033f9734-12e2-465f-873d-19bcd978b464


# ╔═╡ e80b1964-d744-4a1b-b585-0510feadddb4
md""" ## plot """

# ╔═╡ 7a68c8ba-06e5-41f3-b4d1-d9666c7db85d


# ╔═╡ 8927d7d0-0b9e-482d-b816-6492af6902d6
md""" ## save to binary file"""

# ╔═╡ 246ad0d9-1a06-4cad-b33e-141a841574df


# ╔═╡ c3e3bcab-6f5c-466d-82c5-071737630be6
function write_bin(filename::String, A::Array{Float32})
    open(filename, "w") do io
        # hton: convert to big-endian explicitly
        for x in A
            write(io, hton(x))
        end
    end
end

# ╔═╡ 91bd87af-390c-4f37-a8c3-73785fb4d75f
function nearest_grid_index(x, y, XC, YC)
    # circular distance in longitude, handles wraparound at ±180
    dlon = XC .- x
    dlat = YC .- y
    d2 = dlon.^2 .+ dlat.^2
    return argmin(d2)   # returns CartesianIndex{2}
end

# ╔═╡ e31eabe5-d495-4c22-9f7f-b79e594b9093
begin

    XC = write(Γ6.XC)
    YC = write(Γ6.YC)
    # target
    XC_tar = write(Γ.XC)
    YC_tar = write(Γ.YC)

    # for each regional point on ECCO grid, find the nearest point on IAP, then delete the duplicate ones
    region_masks = [0*Γ.XC for i in 1:nregs]
    for ireg = 1:nregs
        # masked points on ECCO grid
        mask_ref = write(region_masks_ref[ireg])
        id = findall(mask_ref .== 1)
        x_ref = XC[id] #
        y_ref = YC[id]
        np_ref = length(id)
        
        # find the nearest points on GLORYS grid
        idxs = Vector{CartesianIndex{2}}(undef, np_ref)
        for i = 1:np_ref
            idxs[i] = nearest_grid_index(x_ref[i], y_ref[i], XC_tar, YC_tar)
        end
        # remove duplicates
        unique_idxs = unique(idxs)
        depth_vals = write(Γ.Depth)[unique_idxs]
        unique_idxs = unique_idxs[depth_vals .!= 0.]
        println("ECCO points: ", np_ref, "  unique grid points on GLORYS: ", 
                length(unique_idxs))
        # ij_pairs = Tuple.(unique_idxs)          # Vector{Tuple{Int,Int}}
        unique_lonlat = [(XC_tar[idx], YC_tar[idx]) for idx in unique_idxs]
        mask2d = zeros(Float32, size(XC_tar))
        mask2d[unique_idxs] .= 1.
        
        # println(findall(write(Γ.Depth)[unique_idxs] .== 0.))
        region_masks[ireg] = γ.read(mask2d,Γ.XC)
    end
end

# ╔═╡ 11bfab01-690a-4dac-820c-af7b01056a57
begin
	fig = Figure()
    # data = write(Γ.Depth)
    data = write(Γ.hFacC)[:,:,1]
    
	ax = Axis(fig[1, 1])
    sc = scatter!(
        ax,
        XC_tar[:],YC_tar[:],
        color = data[:],
        colormap = :berlin,
        markersize = 2,
        # colorrange = (0,6000),
    )

	fig
end

# ╔═╡ 28b5eee9-3966-41b2-b4d4-365ced815242
let
    data = write(Γ.Depth) #* μ 
    reg = 1
    
    fig = Figure()   # size = (600, 500)

    ax = Axis(fig[1, 1],titlesize = 20, title = "mask for GLORYS: $(region_names[reg])")
 #    sc = scatter!(
 #        ax,
 #        XC_tar[:],YC_tar[:],
 #        color = data[:],
 #        colormap = :berlin,
 #        markersize = 12,
 #        colorrange = (0,6000),
 #    )
	# Colorbar(fig[1, 2], sc)

    #---- plot mask
    id = findall(data .== 0) # = unique_idxs
    scatter!(ax,XC_tar[id], YC_tar[id], color=:grey, markersize=2)
    
    id = findall(write(region_masks[reg]) .== 1) # = unique_idxs
    println(minimum(write(Γ.Depth)[id]))
    # println(id)
    scatter!(ax,XC_tar[id], YC_tar[id], color=:red, markersize=4)
    
    xlims!(ax, -85, 5)
    ylims!(ax, 20, 90)	

    fig
end

# ╔═╡ 4d1e8d0e-033d-4036-b5ff-190a216ced92
begin
    to_array(x) = x
    to_array(x::MeshArrays.AbstractMeshArray) = MeshArrays.write(x)
    
    outdir = joinpath(homedir(), "repos/CMEMS_DTOR_Hackathon_2026/data")
    mkpath(outdir)
    
    specs = Tuple("$(name)_maskC_GLORYS" => mask for (name, mask) in zip(region_names, region_masks))
    
    # keep converted arrays for plotting
    mask_arrays = Dict{String,Any}()
    
    for (fname, mask) in specs
        arr = to_array(mask)          # keep regular Array version
        mask_arrays[fname] = arr      # store for plotting
    
        fnm = joinpath(outdir, fname)
        if !isfile(fnm)
            write_bin(fnm, Float32.(arr))
        end
    end
end

# ╔═╡ 6e4c8316-3d6f-4e27-86ed-4ea381bef0a5
begin
    """
    Return grid-cell edges from regularly spaced cell-center coordinates.
    """
    function cell_edges(x::AbstractVector)
        Δx = diff(x)

        @assert length(x) ≥ 2
        @assert all(isapprox.(Δx, first(Δx); atol=1e-4, rtol=1e-4)) ||
            error("grid_factors assumes regularly spaced coordinates.")

        return vcat(
            first(x) - first(Δx) / 2,
            x[1:end-1] .+ first(Δx) / 2,
            last(x) + first(Δx) / 2,
        )
    end

    """
    Create a one-face MeshArrays latitude-longitude grid from
    GLORYS longitude and latitude cell-center vectors.

    Returns:
        Γ  : grid fields, including XC, YC, XG, YG, RAC, dxC, dyC, ...
        γ  : underlying MeshArrays gcmgrid
    """
    function glorys_grid(lon::AbstractVector, lat::AbstractVector)

        # Ensure standard floating-point coordinate vectors.
        xc = Float64.(lon)
        yc = Float64.(lat)

        @assert issorted(xc)
        @assert issorted(yc)

        # Grid-cell-edge coordinates.
        xg = cell_edges(xc)
        yg = cell_edges(yc)

        # This is the input expected by grid_factors / GridSpec_default.
        xy = (; xc, yc, xg, yg)

        # Creates XC, YC, XG, YG, RAC, dxC, dyC, etc.
        Γ = MeshArrays.Grids_simple.grid_factors(xy)

        # The underlying one-face MeshArrays grid specification.
        γ = Γ.XC.grid

        return γ, Γ
    end
end;

# ╔═╡ Cell order:
# ╠═6c84e1e0-a9b1-45a1-b499-5df954d1c6ef
# ╠═41809153-52b3-4f91-8edd-79e5522bad36
# ╠═447ba694-6a19-4400-b0da-79e97b8f5c51
# ╠═f055de04-692e-40d4-8219-3b1d8a95a2cc
# ╠═db7654cb-b219-4216-ae63-cd9900f662a8
# ╠═760143e6-e022-4285-ae1e-2a058cbcc749
# ╠═e546615b-ea5d-4203-964b-de68e2a0d65e
# ╠═921a9e4e-f4b1-430c-9884-5619501a5337
# ╠═1a6f854d-f49e-4603-8652-f467488f5d89
# ╠═dedef362-5070-4901-b322-a3d2894968f7
# ╠═8c4e5825-f5de-4859-9b82-cbdab561d5a5
# ╠═cd56d103-def3-438c-a37f-e28d389df823
# ╠═f4668412-0053-44d4-8ba4-9eed43e7c47a
# ╠═11bfab01-690a-4dac-820c-af7b01056a57
# ╠═3c9add6e-62ed-44f1-b8ff-ea50b9991873
# ╠═e84f674a-4e24-48cd-8dbb-53c345f71ee9
# ╠═5f511f7b-4cff-4691-961b-21b619319a33
# ╠═2f6bd8b9-1123-43a8-a769-80fb2139dbb0
# ╠═033f9734-12e2-465f-873d-19bcd978b464
# ╠═e31eabe5-d495-4c22-9f7f-b79e594b9093
# ╠═e80b1964-d744-4a1b-b585-0510feadddb4
# ╠═28b5eee9-3966-41b2-b4d4-365ced815242
# ╠═7a68c8ba-06e5-41f3-b4d1-d9666c7db85d
# ╠═8927d7d0-0b9e-482d-b816-6492af6902d6
# ╠═4d1e8d0e-033d-4036-b5ff-190a216ced92
# ╠═246ad0d9-1a06-4cad-b33e-141a841574df
# ╠═c3e3bcab-6f5c-466d-82c5-071737630be6
# ╠═91bd87af-390c-4f37-a8c3-73785fb4d75f
# ╠═6e4c8316-3d6f-4e27-86ed-4ea381bef0a5
