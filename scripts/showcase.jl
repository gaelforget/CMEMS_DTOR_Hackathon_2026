### A Pluto.jl notebook ###
# v1.0.1

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 02442c98-2d23-4426-8c9a-ff04cf3f0d17
begin
	import Pkg  
    Pkg.activate(joinpath(@__DIR__, ".."))
end

# ╔═╡ 368c2d94-5bec-4325-a72b-2946269e6cb2
using PlutoUI, JLD2, CairoMakie, Statistics, Dates

# ╔═╡ 9757f64d-2ae0-4be4-8e33-df19a650f189
begin
	using Markdown
	using InteractiveUtils
end

# ╔═╡ 02fe9d6e-d494-11f1-8dca-13f8497eb318
md""" # Labrador Sea stratification viewer

Choose a dataset with the tab (radio buttons), then a variable. 

"""

# ╔═╡ 619bd532-5ee0-404c-a5e5-858a062b3081
md"""
## Choose dataset
$(@bind dataset Radio(["ECCOv4r6" => "ECCOv4r6", "GLORYS12V1" => "GLORYS12V1"], default="ECCOr6"))
"""

# ╔═╡ 5184f12f-3755-44f6-ad91-1a64f2136ecc
md""" 
## Regional mask
"""

# ╔═╡ 42f7d40a-86d1-40f2-8c48-4221d8d3860c
md""" 
## Vertical structures 
"""

# ╔═╡ 966372e8-e00b-41d5-82de-eec662be9efa
function plot_strat_vert(ds, dataset)
	T2d = ds["T_al"]
	S2d = ds["S_al"]
	n2d = ds["N2"]              
	t_diag = ds["t_diag"]     
	zc = ds["zc"]
	logN = map(x -> x > 0 ? log10(x) : NaN, n2d)
	tnum = Dates.value.(t_diag .- t_diag[1]) ./ 86_400_000   

	fig = Figure(size = (600, 800))
	ax  = Axis(fig[1, 1], xlabel = "", ylabel = "Depth [m]",
	           title = "Temperature in Labrador Sea ($dataset)",
          	# fontsize = 16,
	           xticks = (tnum[1:6:end], Dates.format.(t_diag[1:6:end], "yyyy-mm"))
			  )

	hm = heatmap!(ax, tnum, zc, T2d';
				  colormap = :thermal,
				colorrange = (2,10))
	# contours on top of the heatmap
	levels = (4:2:8)
	contour!(ax, tnum, zc, T2d'; levels = levels,
	         color = :black, linewidth = 0.7,
	         labels = true, labelsize = 12)

	Colorbar(fig[1, 2], hm, label = "[°C]")
	ylims!(ax, -2000, 0)
	
	#--------------- middle S
	ax2  = Axis(fig[2, 1], xlabel = "", ylabel = "Depth [m]",
	           title = "Salinity in Labrador Sea ($dataset)",
	           xticks = (tnum[1:6:end], Dates.format.(t_diag[1:6:end], "yyyy-mm"))
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
	           title = "Buoyancy frequency N² in the Labrador Sea ($dataset)",
	           xticks = (tnum[1:6:end], Dates.format.(t_diag[1:6:end], "yyyy-mm"))
			   )
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

# ╔═╡ ea3ddfe9-7c74-4587-9efd-0f9760ad32ff


# ╔═╡ c9b20d80-2a0b-448e-a512-a78daa5977fa
begin
	fnm1  = joinpath(homedir(), "repos/CMEMS_DTOR_Hackathon_2026/data/strat_Labrador_ECCOr6.jld2")
	fnm2 = joinpath(homedir(), "repos/CMEMS_DTOR_Hackathon_2026/data/strat_Labrador_GLORYS.jld2")
	files = Dict("ECCOv4r6" => fnm1, "GLORYS12V1" => fnm2)
end

# ╔═╡ 30ba191c-db10-4bc8-84f7-0626cf2d44ac
begin
	path = files[dataset]
	ds = load(path)            # Dict{String,Any} of everything in the JLD2 file
	varnames = sort(collect(keys(ds)))
	md"Loaded $(path) with $(length(varnames)) entries."
end

# ╔═╡ 212db2e8-7f49-4cf8-8f9b-4f052933274e
begin
	depth = ds["depth2d"]
	XC = ds["XC"]
	YC = ds["YC"]
	mask = ds["mask2d"]
end;

# ╔═╡ c68080b3-c87f-433a-b43f-4c939aa4c21b
let
    fig = Figure()   

    ax = Axis(fig[1, 1],titlesize = 20, title = "Labrador Sea mask, $dataset")

    id = findall(depth .== 0) 
    scatter!(ax,XC[id], YC[id], color=:grey, markersize=4)
    
    id = findall(mask .== 1) 
    scatter!(ax,XC[id], YC[id], color=:red, markersize=4)
    
    xlims!(ax, -85, 5)
    ylims!(ax, 20, 90)	

    fig
end

# ╔═╡ c85c0bc2-2c24-4785-8060-f663f933fa57
fig = plot_strat_vert(ds,dataset);

# ╔═╡ 183f2708-bbcc-4649-bb63-96ea5cc481aa
fig

# ╔═╡ Cell order:
# ╠═02fe9d6e-d494-11f1-8dca-13f8497eb318
# ╠═02442c98-2d23-4426-8c9a-ff04cf3f0d17
# ╠═368c2d94-5bec-4325-a72b-2946269e6cb2
# ╠═9757f64d-2ae0-4be4-8e33-df19a650f189
# ╠═183f2708-bbcc-4649-bb63-96ea5cc481aa
# ╠═619bd532-5ee0-404c-a5e5-858a062b3081
# ╠═30ba191c-db10-4bc8-84f7-0626cf2d44ac
# ╠═5184f12f-3755-44f6-ad91-1a64f2136ecc
# ╠═212db2e8-7f49-4cf8-8f9b-4f052933274e
# ╠═c68080b3-c87f-433a-b43f-4c939aa4c21b
# ╠═42f7d40a-86d1-40f2-8c48-4221d8d3860c
# ╠═c85c0bc2-2c24-4785-8060-f663f933fa57
# ╠═966372e8-e00b-41d5-82de-eec662be9efa
# ╠═ea3ddfe9-7c74-4587-9efd-0f9760ad32ff
# ╠═c9b20d80-2a0b-448e-a512-a78daa5977fa
