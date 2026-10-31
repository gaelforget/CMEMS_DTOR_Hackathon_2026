### A Pluto.jl notebook ###
# v1.0.1

using Markdown
using InteractiveUtils

# ╔═╡ 0bdf8a98-184d-4345-9cc8-b4062d9298eb
begin
	import Pkg  
    Pkg.activate(joinpath(@__DIR__, ".."))
end

# ╔═╡ cb7c3710-4fe9-4652-9f8e-e1b97a81e89e
using MeshArrays, Glob, Statistics, NCDatasets, Dates, CairoMakie, JLD2, MITgcm, GibbsSeaWater

# ╔═╡ d370dc16-d86f-49fc-bbd4-092b973b35b8
using LinearAlgebra: dot

# ╔═╡ 66de41e0-d3e5-11f1-868d-b3b12683fd6e
md""" # read and compute horizontal means of T/S fields from GLORYS"""

# ╔═╡ cf5f04a2-943f-46d4-ae51-7c1c79b2f5a1
begin
	grpdir = get(ENV, "GROUP", "")
	datadir = joinpath(grpdir, "GLORYS/GLOBAL_MULTIYEAR_PHY_001_030/cmems_mod_glo_phy_my_0.083deg_P1M-m_202311")
	println(datadir)
end

# ╔═╡ c1de377f-0601-4d4d-a3a0-409d7495b506
md""" #### get all filenames & read some of them"""

# ╔═╡ fbe38db6-6bfa-4219-ac7b-d83b5154f585
begin
    files = String[]
    
    for year in 2023:2025
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

# ╔═╡ bdebe707-40cc-4f2b-a2af-f8cfdb83e614
ds = NCDataset(files[1:1], "r";
    aggdim = "time",
    deferopen = true,
)

# ╔═╡ 9c88be6f-1f4c-4cc5-8bc4-52d95cea201e
md""" ## Grid - build using MeshArrays """

# ╔═╡ 72f9481a-2e63-476e-b8c4-30b43b45e43c
md""" ## read mask """

# ╔═╡ 510798a2-9f27-4f82-9b27-5293460768be
mask_fnm = joinpath(homedir(), 
					"repos/CMEMS_DTOR_Hackathon_2026/data/lab_maskC_GLORYS")

# ╔═╡ a35fc829-7688-4cf6-9c29-d7b217ae1cf7
md""" ## compute area mean"""

# ╔═╡ 330d3fd2-10d3-4c89-9545-cddf44c086e1
md""" ## compute density and ∂ρ/∂z from using EOS"""

# ╔═╡ d5c6ca2f-a513-4c6b-aaa3-272e2f4b070e
md""" ## plot """

# ╔═╡ 107c3b37-3cbb-4138-9501-07c664f0bd3c
md""" ## save to JLD2 """

# ╔═╡ be9bf5e6-a294-49a5-9587-bb6cf6af3aa2
savefnm = joinpath(homedir(), 
				   "repos/CMEMS_DTOR_Hackathon_2026/data/strat_Labrador_GLORYS.jld2")

# ╔═╡ 2b65aa9c-6b8b-4d02-b06c-9190265aa558


# ╔═╡ 5dae2e8f-b316-4050-b5cd-99cd79bf5aa2
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

# ╔═╡ 88395a56-65c0-4adc-bcd2-b60943c3f393
begin
	RC = ds["depth"][:] # depth of layer center
	lon = ds["longitude"][:]
	lat = ds["latitude"][:]
	temp = ds["thetao"][:,:,:,1]
	γ, Γ_raw = glorys_grid(lon, lat)
	close(ds)
end

# ╔═╡ 8f68592f-ba06-4f5b-b7e5-c2e038178b59
begin
	RF = [0; 0.5*(RC[1:end-1]+RC[2:end]); RC[end]+0.5*(RC[end]-RC[end-1])] 
	DRF = diff(RF)
	hFacC = read(1.0*(!ismissing).(temp),γ)
	depth = Γ_raw.XC * 0
	for t in axes(hFacC, 1), k in axes(hFacC, 2)
		depth[t] += hFacC[t, k] .* DRF[k]
	end
	Γ = merge(Γ_raw, (RC = RC, RF = RF, DRF = DRF, hFacC = hFacC, Depth=depth))
	nx, ny, nz = size(write(Γ.hFacC))
end

# ╔═╡ 1c842f6d-0dbf-4404-b175-895471d7c830
read_mask(fname) = begin
    raw = Array{Float32}(undef, nx, ny)
    MITgcm.read!(fname, raw)
    raw .= ntoh.(raw)
end

# ╔═╡ ef7a154e-58e8-499f-96ec-c43ebf68ba36
mask2d = read_mask(mask_fnm)

# ╔═╡ 85b9362e-a149-4328-ae35-1fdd6684c573
begin
    XC = write(Γ.XC)
    YC = write(Γ.YC)
    depth2d = write(Γ.Depth)
end

# ╔═╡ d98f5ace-96e7-4e79-ae8e-ddcd5b625996
let
    data = write(Γ.Depth) #* μ 
    reg = 1

    fig = Figure()   # size = (600, 500)

    ax = Axis(fig[1, 1],titlesize = 20, title = "mask for GLORYS")

    id = findall(data .== 0) # = unique_idxs
    scatter!(ax,XC[id], YC[id], color=:grey, markersize=2)
    
    id = findall(mask2d .== 1) # = unique_idxs
    scatter!(ax,XC[id], YC[id], color=:red, markersize=4)
    
    xlims!(ax, -85, 5)
    ylims!(ax, 20, 90)	

    fig
end

# ╔═╡ 321e6703-b30f-4b6e-9468-c467692ca126
begin
    # area weights inside the region (plain Float64 matrix, nx × ny)
    w2d = Float64.(mask2d) .* Float64.(write(Γ.RAC))

    # bounding box of the region, so we only read that part of each file
    ci = findall(>(0), w2d)
    i1, i2 = extrema(c[1] for c in ci)
    j1, j2 = extrema(c[2] for c in ci)
    w_sub = w2d[i1:i2, j1:j2]
end

# ╔═╡ 40f1ff1a-5f73-4cc0-bc14-ca56f543693e
zc =  -Γ.RC

# ╔═╡ 75c95477-df92-4710-8522-ac184cd2a51d
# area-weighted mean per layer; skips missing/NaN points (land, below seafloor)
function layer_mean!(out, F, w)
    nx, ny, nz = size(F)
    @inbounds for k in 1:nz
        num = 0.0
        den = 0.0
        for j in 1:ny, i in 1:nx
            wij = w[i, j]
            wij == 0 && continue
            x = F[i, j, k]
            if !ismissing(x) && !isnan(x)
                num += x * wij
                den += wij
            end
        end
        out[k] = den > 0 ? num / den : NaN
    end
    return out
end

# ╔═╡ 81177436-37bb-4178-8494-7abb0ec0c39d
function process_files(files, it_do, nz, i1, i2, j1, j2, w_sub)
    nt_do  = length(it_do)
    T_al   = Matrix{Float64}(undef, nz, nt_do)
    S_al   = Matrix{Float64}(undef, nz, nt_do)
    t_diag = Vector{DateTime}(undef, nt_do)

    for (it, f) in enumerate(files[it_do])
        # @info "processing file $it / $nt_do" _id=:fileprog
        NCDataset(f, "r") do ds
            t_diag[it] = ds["time"][1]
            T3d = ds["thetao"][i1:i2, j1:j2, :, 1]   # keeps missing, no copy
            S3d = ds["so"][i1:i2, j1:j2, :, 1]
            layer_mean!(view(T_al, :, it), T3d, w_sub)
            layer_mean!(view(S_al, :, it), S3d, w_sub)
        end                                          # file closed here
    end
    return T_al, S_al, t_diag
end

# ╔═╡ 7e746f9b-6780-48e9-a978-c0169f073951
begin
    it_do = 1:length(files)          # use 1:1 for a quick test
    T_al, S_al, t_diag = @time process_files(files, it_do, nz, i1, i2, j1, j2, w_sub)
end

# ╔═╡ 9e7f7252-6e47-40db-a0cb-87c75ccc7121
function ddz_center(f, zc)
    nz, nt = size(f)
    out = fill(NaN, nz, nt)
    @inbounds for k in 2:nz-1
        h1 = zc[k]   - zc[k-1]
        h2 = zc[k+1] - zc[k]
        a  = -h2 / (h1 * (h1 + h2))
        b  = (h2 - h1) / (h1 * h2)
        c  =  h1 / (h2 * (h1 + h2))
        for t in 1:nt
            out[k, t] = a * f[k-1, t] + b * f[k, t] + c * f[k+1, t]
        end
    end
    # one-sided at the top and bottom
    out[1,  :] .= (f[2, :]  .- f[1, :])    ./ (zc[2]  - zc[1])
    out[nz, :] .= (f[nz, :] .- f[nz-1, :]) ./ (zc[nz] - zc[nz-1])
    return out
end

# ╔═╡ f5365f03-5e50-4287-9a35-ed2b01378f64
begin
	g, ρ0 = 9.81, 1025.0
	rho_al = similar(T_al)
	dσdz = similar(T_al)
	N2 = similar(T_al)
	lon2d = 310 .* ones(size(T_al))
	lat2d = 60 .* ones(size(T_al))
	z2d = repeat(-Γ.RC, 1, size(T_al,2)) 
	S2d = S_al
	T2d = T_al
	p2d = gsw_p_from_z.(z2d, lat2d)
	SA = gsw_sa_from_sp.(S2d, p2d, lon2d, lat2d)
	CT = gsw_ct_from_pt.(SA, T2d)
	rho_al = gsw_sigma0.(SA, CT)
	dσdz = ddz_center(rho_al, -Γ.RC)

	N2 = -(g / ρ0) .* dσdz
end

# ╔═╡ ba203d3c-a7b3-401d-a862-8f39bca1cb6a
let
	f2d = T_al
	S2d = S_al
	n2d = N2              
	logN = map(x -> x > 0 ? log10(x) : NaN, n2d)
	
	tnum = Dates.value.(t_diag .- t_diag[1]) ./ 86_400_000   # days since the start
	# tnum = it_do
		
	fig = Figure(size = (800, 1000))
	ax  = Axis(fig[1, 1], xlabel = "", ylabel = "Depth [m]",
	           title = "Area-mean temperature in Labrador Sea (GLORYS)",
	           xticks = (tnum[1:4:end], Dates.format.(t_diag[1:4:end], "yyyy-mm"))
			  )

	hm = heatmap!(ax, tnum, zc, f2d';
				  colormap = :thermal,
				colorrange = (2,10))
	# contours on top of the heatmap
	# levels = floor(minimum(filter(!isnan, f2d))):1:ceil(maximum(filter(!isnan, f2d)))
	levels = (4:2:8)
	contour!(ax, tnum, zc, f2d'; levels = levels,
	         color = :black, linewidth = 0.7,
	         labels = true, labelsize = 12)

	Colorbar(fig[1, 2], hm, label = "[°C]")
	ylims!(ax, -2000, 0)
	
	#--------------- middle S
	ax2  = Axis(fig[2, 1], xlabel = "Time", ylabel = "Depth [m]",
	           title = "Area-mean salinity in Labrador Sea (GLORYS)",
	           xticks = (tnum[1:4:end], Dates.format.(t_diag[1:4:end], "yyyy-mm"))
			   )
	hm2 = heatmap!(ax2, tnum, zc, S2d';
				  colormap = (:haline),
				colorrange = (34.3,34.9)
				 )
	levels = [34.6, 34.8, 34.9]
	contour!(ax2, tnum, zc, S2d'; levels = levels,
	         color = :black, linewidth = 0.7,
	         labels = true, labelsize = 12)
	Colorbar(fig[2, 2], hm2, label = "[PSU]")
	ylims!(ax2, -2000, 0)

	#--------------- botom N^2
	ax3  = Axis(fig[3, 1], xlabel = "", ylabel = "Depth [m]",
	           title = "Buoyancy frequency N² in the Labrador Sea (GLORYS)",
	           xticks = (tnum[1:4:end], Dates.format.(t_diag[1:4:end], "yyyy-mm")))
	hm = heatmap!(ax3, tnum, zc, logN';
	              colormap = :deep,
	              colorrange = (-8, -3),           # 1e-5 to 1e-2 kg m⁻⁴
	              lowclip = :white, highclip = :black,
	              nan_color = :lightgray)           
	contour!(ax3, tnum, zc, logN'; levels = [-6, -4],
	         color = :black, linewidth = 0.8, labels = true, labelsize = 12)
	Colorbar(fig[3, 2], hm, label = "log₁₀(N²)  [s⁻²]",
	         ticks = -8:-3)
	ylims!(ax3, -2000, 0)
	
	fig
end

# ╔═╡ 3c17585f-54b6-4c30-9e71-211ce90634b6
if !isfile(savefnm)
    jldsave(savefnm;
        zc, t_diag, T_al, S_al, N2, XC, YC, depth2d, mask2d)
end

# ╔═╡ Cell order:
# ╠═66de41e0-d3e5-11f1-868d-b3b12683fd6e
# ╠═0bdf8a98-184d-4345-9cc8-b4062d9298eb
# ╠═cb7c3710-4fe9-4652-9f8e-e1b97a81e89e
# ╠═d370dc16-d86f-49fc-bbd4-092b973b35b8
# ╠═cf5f04a2-943f-46d4-ae51-7c1c79b2f5a1
# ╠═c1de377f-0601-4d4d-a3a0-409d7495b506
# ╠═fbe38db6-6bfa-4219-ac7b-d83b5154f585
# ╠═bdebe707-40cc-4f2b-a2af-f8cfdb83e614
# ╠═9c88be6f-1f4c-4cc5-8bc4-52d95cea201e
# ╠═88395a56-65c0-4adc-bcd2-b60943c3f393
# ╠═8f68592f-ba06-4f5b-b7e5-c2e038178b59
# ╠═72f9481a-2e63-476e-b8c4-30b43b45e43c
# ╠═510798a2-9f27-4f82-9b27-5293460768be
# ╠═1c842f6d-0dbf-4404-b175-895471d7c830
# ╠═ef7a154e-58e8-499f-96ec-c43ebf68ba36
# ╠═85b9362e-a149-4328-ae35-1fdd6684c573
# ╠═d98f5ace-96e7-4e79-ae8e-ddcd5b625996
# ╠═a35fc829-7688-4cf6-9c29-d7b217ae1cf7
# ╠═321e6703-b30f-4b6e-9468-c467692ca126
# ╠═7e746f9b-6780-48e9-a978-c0169f073951
# ╠═330d3fd2-10d3-4c89-9545-cddf44c086e1
# ╠═f5365f03-5e50-4287-9a35-ed2b01378f64
# ╠═40f1ff1a-5f73-4cc0-bc14-ca56f543693e
# ╠═d5c6ca2f-a513-4c6b-aaa3-272e2f4b070e
# ╠═ba203d3c-a7b3-401d-a862-8f39bca1cb6a
# ╠═107c3b37-3cbb-4138-9501-07c664f0bd3c
# ╠═be9bf5e6-a294-49a5-9587-bb6cf6af3aa2
# ╠═3c17585f-54b6-4c30-9e71-211ce90634b6
# ╠═2b65aa9c-6b8b-4d02-b06c-9190265aa558
# ╠═5dae2e8f-b316-4050-b5cd-99cd79bf5aa2
# ╠═75c95477-df92-4710-8522-ac184cd2a51d
# ╠═81177436-37bb-4178-8494-7abb0ec0c39d
# ╠═9e7f7252-6e47-40db-a0cb-87c75ccc7121
