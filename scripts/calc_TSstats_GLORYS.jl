### A Pluto.jl notebook ###
# v1.0.1

using Markdown
using InteractiveUtils

# ╔═╡ 66b9a22e-5c7b-4bba-8dd7-4ad56384e7ee
begin
	import Pkg  
    Pkg.activate(joinpath(@__DIR__, ".."))
end

# ╔═╡ 9a275f2d-6add-4053-bada-51bfa8e005ce
using MeshArrays, Glob, Statistics, NCDatasets, Dates, CairoMakie, JLD2

# ╔═╡ 93515a8a-ccbf-11f1-b95c-85ef83f600db
md""" # read and compute zonal means of T/S fields from GLORYS"""

# ╔═╡ 30b7c70b-8039-4e94-91ac-b2f7ea7fd472
grpdir = get(ENV, "GROUP", "")

# ╔═╡ 8af64a79-31e9-49af-af13-dcf9c26cb2e9
begin
	datadir = joinpath(grpdir, "GLORYS/GLOBAL_MULTIYEAR_PHY_001_030/cmems_mod_glo_phy_my_0.083deg_P1M-m_202311")
	println(datadir)
end

# ╔═╡ a0a96a34-3c05-4ad2-8891-0be2919a16a0
md""" #### get all filenames"""

# ╔═╡ 5ce85ab9-a2be-4314-932b-26fa00f3e187
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
    println.(files)
end

# ╔═╡ a75b3ba2-9110-420e-b71d-fa88452416aa
md""" #### read data using NCDataset vs YAXArrays"""

# ╔═╡ c6b76932-a0d2-4b8f-bb36-e0e7cf756aa3
ds = NCDataset(files, "r";
    aggdim = "time",
    deferopen = true,
)

# ╔═╡ aa8669cd-6477-4f28-8d62-53d2c680a82d
# import NetCDF

# ╔═╡ bab62b50-66df-45c8-8bc7-37689249775d
# ds_yax = open_dataset(files[1], driver=:netcdf)

# ╔═╡ 3edc5492-3a3d-4641-89b5-09816d88549c
md""" ## Grid - build using MeshArrays """

# ╔═╡ 37a026db-6c70-4abc-9d9a-6a42255474e4
begin
	lon = ds["longitude"][:]
	lat = ds["latitude"][:]
	depth = ds["depth"][:]
	γ, Γ = glorys_grid(lon, lat)
	# γ4 = GridSpec(ID=:IAP)
	# Γ4 = GridLoad(γ4;option="full")
end

# ╔═╡ ad502098-1e42-4760-9cf4-cbd2a1ad5dec
begin
	println(extrema(Γ.XC[1,1][:,222] .- lon))
	println(extrema(Γ.YC[1,1][333,:] .- lat))
	size(Γ.XC[1])
end

# ╔═╡ 401bcd6c-436c-493f-81aa-ebb88ecec00a
md""" ## Zonal mean """

# ╔═╡ 05dfa1d6-7094-4238-b259-fc2924ccc01d
begin
	v = ds["thetao"]   
	nx, ny, nz, nt = size(v)
	ik = 1:nz
	it = [1]
	nz_do = length(ik)
	nt_do = length(it)
	
	times = ds["time"][:]         
	ti = times[it]
	# close(ds)
end

# ╔═╡ d6d4fd2a-00ab-4bef-ace7-183ca43d6b2a
@time begin
    t_zm = Array{Float64}(undef, ny, nz_do, nt_do)  
    for t in it
        slice = v[:, :, :, t]                         
        t_zm[:, :, t] = nanzonalmean(coalesce.(slice, NaN))
        
        if t % 20 == 0 || t == nt
            println("Step $t / $nt")
        end
    end
end

# ╔═╡ 0da6fd03-2523-4bc8-99dd-02f211830b0e
md""" ### plot lat-depth of T_zm at one time """

# ╔═╡ 7b712a9a-da04-49da-affc-2e79bb51c7d9
let
    it = 1
    fplt = t_zm[:, :, it]
    
    vals = filter(!isnan, fplt)
    println(extrema(vals))
    month_label = Dates.format(ti[it], "yyyy-mm")

    fig, ax, hm = heatmap(lat, depth, fplt;
            colormap = :thermal,
            colorrange = (0, 30),       
            axis = (xlabel="Latitude", ylabel="Depth (m)", yreversed=true, title="Zonal mean θ at $(month_label)"))
    Colorbar(fig[1, 2], hm, label="°C")
    # ylims!(4000,0)
    fig
end

# ╔═╡ 33ee3171-8922-4b68-87df-7cc6cf069d39
md""" ## *Functions* """

# ╔═╡ 7ce0b78b-9460-4de8-9a14-f6cb95fbe581
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

# ╔═╡ 090165a6-714f-496f-8603-f3c8a8267c80


# ╔═╡ 9a49cda2-95c0-4db2-915b-a475657e0068
function nanzonalmean(A)           # A: (lon, lat, depth)
    valid = .!isnan.(A)
    cnt   = sum(valid; dims=1)
    tot   = sum(ifelse.(valid, A, 0f0); dims=1)
    return dropdims(tot ./ cnt; dims=1)   # 0/0 gives NaN over land
end

# ╔═╡ Cell order:
# ╠═93515a8a-ccbf-11f1-b95c-85ef83f600db
# ╠═66b9a22e-5c7b-4bba-8dd7-4ad56384e7ee
# ╠═9a275f2d-6add-4053-bada-51bfa8e005ce
# ╠═30b7c70b-8039-4e94-91ac-b2f7ea7fd472
# ╠═8af64a79-31e9-49af-af13-dcf9c26cb2e9
# ╠═a0a96a34-3c05-4ad2-8891-0be2919a16a0
# ╠═5ce85ab9-a2be-4314-932b-26fa00f3e187
# ╠═a75b3ba2-9110-420e-b71d-fa88452416aa
# ╠═c6b76932-a0d2-4b8f-bb36-e0e7cf756aa3
# ╠═aa8669cd-6477-4f28-8d62-53d2c680a82d
# ╠═bab62b50-66df-45c8-8bc7-37689249775d
# ╠═3edc5492-3a3d-4641-89b5-09816d88549c
# ╠═37a026db-6c70-4abc-9d9a-6a42255474e4
# ╠═ad502098-1e42-4760-9cf4-cbd2a1ad5dec
# ╠═401bcd6c-436c-493f-81aa-ebb88ecec00a
# ╠═05dfa1d6-7094-4238-b259-fc2924ccc01d
# ╠═d6d4fd2a-00ab-4bef-ace7-183ca43d6b2a
# ╠═0da6fd03-2523-4bc8-99dd-02f211830b0e
# ╠═7b712a9a-da04-49da-affc-2e79bb51c7d9
# ╠═33ee3171-8922-4b68-87df-7cc6cf069d39
# ╠═7ce0b78b-9460-4de8-9a14-f6cb95fbe581
# ╠═090165a6-714f-496f-8603-f3c8a8267c80
# ╠═9a49cda2-95c0-4db2-915b-a475657e0068
