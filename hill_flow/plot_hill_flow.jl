using CairoMakie
using NCDatasets

# =============================================================================
# Animate xz vorticity from hill_flow.jl output.
# Run hill_flow.jl first to produce hill_flow.nc.
# =============================================================================

# --- Seamount geometry (must match hill_flow.jl) ---
const Lx = 20.0
const H  = 2.0
const x₀ = Lx / 3
const h₀ = 0.6H
const σ  = Lx / 10
hill(x) = h₀ * exp(-((x - x₀) / σ)^2) - H

# --- Load data ---
ds = NCDataset("hill_flow.nc")

x = ds["xC"][:]   # cell-center x
z = ds["zC"][:]   # cell-center z
t = ds["time"][:]

# Fields are (x, y=1, z, time); squeeze the Flat y dimension.
ω = dropdims(ds["ω"][:, 1, :, :], dims=2)   # (Nx, Nz, Nt)
u = dropdims(ds["u"][:, 1, :, :], dims=2)

close(ds)

Nt = length(t)

# Mask immersed (below hill) cells with NaN
for k in axes(ω, 2), i in axes(ω, 1)
    if z[k] < hill(x[i])
        ω[i, k, :] .= NaN
        u[i, k, :] .= NaN
    end
end

ω_lim = quantile(filter(!isnan, abs.(vec(ω))), 0.98)

# --- Figure ---
fig = Figure(size=(900, 420))

ax = Axis(fig[1, 1];
          title   = "Vorticity  ω = ∂ᵤu − ∂ₓw",
          xlabel  = "x (m)",
          ylabel  = "z (m)",
          aspect  = DataAspect())

n = Observable(1)
ω_plt = @lift ω[:, :, $n]
title_str = @lift @sprintf("t = %.1f s", t[$n])
Label(fig[0, 1], title_str, fontsize=18)

hm = heatmap!(ax, x, z, ω_plt; colormap=:vik, colorrange=(-ω_lim, ω_lim))
Colorbar(fig[1, 2], hm; label="ω (s⁻¹)")

# Overlay hill profile
x_fine = range(0, Lx, length=500)
z_bottom = hill.(x_fine)
band!(ax, x_fine, fill(-H, length(x_fine)), z_bottom; color=(:gray, 0.8))

# --- Record animation ---
record(fig, "hill_flow.mp4", 1:Nt; framerate=20) do nn
    n[] = nn
end

@info "Animation saved to hill_flow.mp4"
