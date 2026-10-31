### A Pluto.jl notebook ###
# v1.0.1

using Markdown
using InteractiveUtils

# ╔═╡ 5161b8a8-a817-465f-b30c-69437dfd9889
begin
	import Pkg  
    Pkg.activate(joinpath(@__DIR__, ".."))
end

# ╔═╡ b941d29d-c53f-4606-9611-40021699c683
using MeshArrays, Glob, Statistics, NCDatasets, Dates, CairoMakie, JLD2, MITgcm, GibbsSeaWater

# ╔═╡ f53ecfcd-cddc-48e5-b0fc-662a5a1d4ca9
using LinearAlgebra: dot

# ╔═╡ 6372d44a-d3b2-11f1-a733-d96c84d88243
md""" # read and compute horizontal means of T/S fields from ECCOv4r6"""

# ╔═╡ cddb05cd-ae9c-4e8b-980c-4cc7db650850
grpdir = get(ENV, "GROUP", "")

# ╔═╡ 72d34bf0-4f55-45a8-8734-dcb0b2561fb4
begin
	v6dir = "/nfs/micklab005/gforget/2026_datasets/ECCO4r6_output/"
	griddir = joinpath(v6dir,"grid")
	rundir = joinpath(v6dir,"diags_monthly")
	println(rundir)
end

# ╔═╡ 17370f6d-b12c-4f3c-90ad-36fd2e511189
md""" ## Grid """

# ╔═╡ 3ebfe573-35bc-4a57-b48e-ee737bcb2733
begin
	γ = GridSpec("LatLonCap", griddir, ioPrec=Float32)
	Γ = GridLoad(γ;option="full")

	nx, ny, nz = size(write(Γ.hFacC))

	# 1 for oce, NaN for land
	ocean_mask = γ.read(ifelse.(write(Γ.hFacC) .> 0, 1.0, NaN) , Γ.hFacC)
	ocean_mask01 = γ.read(Float64.(write(Γ.hFacC) .> 0), Γ.hFacC)
	# surf
	μ = ocean_mask[:,1] # save as the func: land_mask(Γ)
	μ01 = ocean_mask01[:,1]
	
	zb = -Γ.RF[2:end]
	zc =  Γ.RC
end;

# ╔═╡ 0d99d42b-4d0f-4894-a389-7b12590c21cf
md""" ## 2d mask & 3d volme """

# ╔═╡ ae2b07c3-d18e-4409-b6fd-5780a9f13a3a


# ╔═╡ c63eaaff-1d45-4073-8aee-c5fd6bd4bc66
begin
    # Define region specifications
    mask_fnm = joinpath(expanduser("~/data_small"), "refine_lab_maskC")
    
    # Read all masks efficiently
    read_mask(fname) = begin
        raw = Array{Float32}(undef, nx, ny)
        MITgcm.read!(fname, raw)
        raw .= ntoh.(raw)
        γ.read(raw, Γ.XC)
    end

    # refine_inter_maskC, etc
    region_masks = read_mask(mask_fnm)
    
    # Compute 2d masks
    area2d_masks = similar(Γ.XC)
    
    @inbounds for t in axes(Γ.hFacC, 1)
        v = Γ.RAC[t] 
            area2d_masks[t] = v .* region_masks[t]
    end
    mask2d = write(region_masks)
end

# ╔═╡ f46e832b-b566-47c5-8698-45c962be4436
let
    fig = Figure()    
    data =  Γ.Depth #area2d_masks[1][:] *

    ax = Axis(fig[1, 1],
    )

    # sc = scatter!(
    #     ax,
    #     Γ.XC[:],Γ.YC[:],
    #     color = data[:], # data[:],
    #     # colormap = :balance,
    #     markersize = 6,
    #     # colorrange = (-m, m),
    # )
    
    #---- plot land
    id = findall(isnan.(μ))
    scatter!(ax,Γ.XC[id], Γ.YC[id], color = :grey, markersize=6)

    id = findall((area2d_masks[:]).>0)
    scatter!(ax,Γ.XC[id], Γ.YC[id], color = :red, markersize=6)

    xlims!(ax, -85, 5)
    ylims!(ax, 20, 90)	
    
    fig
end

# ╔═╡ ffdd9578-7cda-48ac-a628-b2fd3bf17350


# ╔═╡ 24beaf66-0807-43fe-a651-00427a974b4d


# ╔═╡ c9fa7731-42f4-4902-a0cb-12d124f19e74
md""" ## read T/S at all times and compute the area means """

# ╔═╡ faed54c1-ab8a-49ca-be75-2a89646c9491
begin
	files_T = glob("THETA_mon_mean.*.data", joinpath(rundir,"THETA_mon_mean"))
	files_S = glob("SALT_mon_mean.*.data", joinpath(rundir,"SALT_mon_mean"))
	@assert length(files_T) == length(files_S)
    nfiles = length(files_T)

    deltaT = 3600.0
    startDate = DateTime(1992, 1, 1, 12, 0, 0)
	it_do = 373:nfiles
	nt_do = length(it_do)
	files_T_do = files_T[it_do]
	files_S_do = files_S[it_do]
end

# ╔═╡ 8ffd26d5-fd2b-41bb-b900-7247d12dd204
md""" ## compute density and ∂ρ/∂z from using EOS"""

# ╔═╡ bc6a94c4-fa9b-46cf-ada1-dd3832418f04


# ╔═╡ 1f5e1fa0-0a6d-4947-9567-02760e8bd8a2
md""" ## plot """

# ╔═╡ 800a8a28-cf0c-453f-81c8-a6ca1393a404
# begin
# 	reg = 1
# 	f2d = T_al[reg]
# 	tnum = Dates.value.(t_diag .- t_diag[1]) ./ 86_400_000   # days since the start

# 	fig = Figure(size = (800, 450))
# 	ax  = Axis(fig[1, 1], xlabel = "Time", ylabel = "Depth [m]",
# 	           title = "Area-mean temperature in Labrador Sea (ECCOv4r6)",
# 			   xticks = (tnum[1:4:end], Dates.format.(t_diag[1:4:end], "yyyy-mm"))
# 			  )
	
# 	hm = heatmap!(ax, tnum, -zc, f2d'; 
# 	              colormap = :thermal)
# 	Colorbar(fig[1, 2], hm, label = "T [°C]")

# 	ylims!(-1000,0)
	
# 	fig
# end

# ╔═╡ 56621b43-7b57-4fa0-a6f3-6a50759c7e5a
md""" ## save to JLD2 """

# ╔═╡ 75bf1bcd-1dd5-4a46-959d-bd04eed2c0f9
savefnm = joinpath(homedir(), 
				   "repos/CMEMS_DTOR_Hackathon_2026/data/strat_Labrador_ECCOr6.jld2")

# ╔═╡ a54d28cb-5ca9-4fb8-a444-f6f1f8a1cb67
begin
	XC = write(Γ.XC)
	YC = write(Γ.YC)
	depth2d = write(Γ.Depth) 
end;

# ╔═╡ 25e96a73-9ff2-4253-88c9-44d5a6f0be38


# ╔═╡ eb20b058-be0f-4a1a-9b4e-cbed3546cc98
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

# ╔═╡ 276a693f-6be7-4b08-9f20-a420e9b56b8c
function compute_T_S_xymean(files_T, files_S, masks, startDate, deltaT, nz;
                                  snapshot_it::Int = length(files_T))
    nfiles = length(files_T)

    t_diag = Vector{DateTime}(undef, nfiles)
    T_al   = Matrix{Float64}(undef, nz, nfiles)
    S_al   = Matrix{Float64}(undef, nz, nfiles)

    maskvecs = vec(write(masks)) 

    T3d_snapshot = nothing
    S3d_snapshot = nothing

    for (it, (fT, fS)) in enumerate(zip(files_T, files_S))
        metaT = read_meta(replace(fT, r"\.data$" => ".meta"))
        metaS = read_meta(replace(fS, r"\.data$" => ".meta"))
        @assert metaT[5] == metaS[5] "T/S timestep mismatch at it=$it"

        t_diag[it] = startDate + Second(round(Int, (metaT[5] - 1) * deltaT))

        T3d = read_mdsio(fT, :THETA)
        S3d = read_mdsio(fS, :SALT)
        for k = 1:nz
            Tvec = vec(T3d[:,:,k])
            Svec = vec(S3d[:,:,k])
            T_al[k,it] = dot(Tvec, maskvecs) / sum(maskvecs)
            S_al[k,it] = dot(Svec, maskvecs) / sum(maskvecs)
        end

        if it == snapshot_it
            T3d_snapshot = T3d   # keep the full 3D field for later analysis
            S3d_snapshot = S3d
        end
    end

    t_mon = Dates.format.(t_diag, dateformat"yyyy-mm")
    return t_diag, t_mon, T_al, S_al, T3d_snapshot, S3d_snapshot
end

# ╔═╡ 19bff648-84eb-4228-a15f-f9dfdc9ed910
@time t_diag, t_mon, T_al, S_al, T3d, S3d =
    compute_T_S_xymean(files_T_do, files_S_do, area2d_masks, startDate, deltaT, nz)

# ╔═╡ c588e995-66b5-427f-b70b-2a2d801beece
begin
	g, ρ0 = 9.81, 1025.0
	rho_al = similar(T_al)
	dσdz = similar(T_al)
	N2 = similar(T_al)
	lon2d = 310 .* ones(size(T_al))
	lat2d = 60 .* ones(size(T_al))
	z2d = repeat(Γ.RC, 1, size(T_al,2)) 
	for id in eachindex(T_al)
		S2d = S_al
		T2d = T_al
	    p2d = gsw_p_from_z.(z2d, lat2d)
		SA = gsw_sa_from_sp.(S2d, p2d, lon2d, lat2d)
		CT = gsw_ct_from_pt.(SA, T2d)
		rho_al = gsw_sigma0.(SA, CT)
		dσdz = ddz_center(rho_al, Γ.RC)

		N2 = -(g / ρ0) .* dσdz
	end
end

# ╔═╡ 12bd81d4-1271-4e1f-8218-f0a0623c03d1
let
	f2d = T_al
	S2d = S_al
	n2d = N2            
	logN = map(x -> x > 0 ? log10(x) : NaN, n2d)

	tnum = Dates.value.(t_diag .- t_diag[1]) ./ 86_400_000   # days since the start

	fig = Figure(size = (800, 1000))
	ax  = Axis(fig[1, 1], xlabel = "", ylabel = "Depth [m]",
	           title = "Area-mean temperature in Labrador Sea (ECCOv4r6)",
	           xticks = (tnum[1:4:end], Dates.format.(t_diag[1:4:end], "yyyy-mm")))

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
	ax2  = Axis(fig[2, 1], xlabel = "", ylabel = "Depth [m]",
	           title = "Area-mean salinity in Labrador Sea (ECCOv4r6)",
	           xticks = (tnum[1:4:end], Dates.format.(t_diag[1:4:end], "yyyy-mm")))
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
	           title = "Buoyancy frequency N² in the Labrador Sea (ECCOv4r6)",
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

# ╔═╡ 9a43ada1-418e-4b85-a28c-b81e6d8dddba
if !isfile(savefnm)
    jldsave(savefnm;
        zc, t_diag, T_al, S_al, N2, XC, YC, depth2d, mask2d) 
end

# ╔═╡ Cell order:
# ╠═6372d44a-d3b2-11f1-a733-d96c84d88243
# ╠═5161b8a8-a817-465f-b30c-69437dfd9889
# ╠═b941d29d-c53f-4606-9611-40021699c683
# ╠═f53ecfcd-cddc-48e5-b0fc-662a5a1d4ca9
# ╠═cddb05cd-ae9c-4e8b-980c-4cc7db650850
# ╠═72d34bf0-4f55-45a8-8734-dcb0b2561fb4
# ╠═17370f6d-b12c-4f3c-90ad-36fd2e511189
# ╠═3ebfe573-35bc-4a57-b48e-ee737bcb2733
# ╠═0d99d42b-4d0f-4894-a389-7b12590c21cf
# ╠═ae2b07c3-d18e-4409-b6fd-5780a9f13a3a
# ╠═c63eaaff-1d45-4073-8aee-c5fd6bd4bc66
# ╠═f46e832b-b566-47c5-8698-45c962be4436
# ╠═ffdd9578-7cda-48ac-a628-b2fd3bf17350
# ╠═24beaf66-0807-43fe-a651-00427a974b4d
# ╠═c9fa7731-42f4-4902-a0cb-12d124f19e74
# ╠═faed54c1-ab8a-49ca-be75-2a89646c9491
# ╠═19bff648-84eb-4228-a15f-f9dfdc9ed910
# ╠═8ffd26d5-fd2b-41bb-b900-7247d12dd204
# ╠═c588e995-66b5-427f-b70b-2a2d801beece
# ╠═bc6a94c4-fa9b-46cf-ada1-dd3832418f04
# ╠═1f5e1fa0-0a6d-4947-9567-02760e8bd8a2
# ╠═800a8a28-cf0c-453f-81c8-a6ca1393a404
# ╠═12bd81d4-1271-4e1f-8218-f0a0623c03d1
# ╠═56621b43-7b57-4fa0-a6f3-6a50759c7e5a
# ╠═75bf1bcd-1dd5-4a46-959d-bd04eed2c0f9
# ╠═a54d28cb-5ca9-4fb8-a444-f6f1f8a1cb67
# ╠═9a43ada1-418e-4b85-a28c-b81e6d8dddba
# ╠═25e96a73-9ff2-4253-88c9-44d5a6f0be38
# ╠═eb20b058-be0f-4a1a-9b4e-cbed3546cc98
# ╠═276a693f-6be7-4b08-9f20-a420e9b56b8c
