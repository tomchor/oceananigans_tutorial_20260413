using Oceananigans
import NCDatasets

# =============================================================================
# Animate hill_flow.jl output.
# Run hill_flow.jl first to produce hill_flow.nc.
# =============================================================================

ωy_ts = FieldTimeSeries("hill_flow.nc", "ωy")
ωz_ts = FieldTimeSeries("hill_flow.nc", "ωz")
w_ts  = FieldTimeSeries("hill_flow.nc", "w")

times = ωy_ts.times
Nt    = length(times)

Ny = size(ωy_ts, 2)   # ωy is Center in y → correct count for xz slice index
Nz = size(ωz_ts, 3)   # ωz is Center in z → correct count for surface slice index

using Statistics: quantile
ωy_lim = quantile(abs.(vec(interior(ωy_ts, :, Ny÷2, :, :))), 0.98)
ωz_lim = quantile(abs.(vec(interior(ωz_ts, :, :, Nz,   :))), 0.98)
w_lim  = quantile(abs.(vec(interior(w_ts,  :, Ny÷2, :, :))), 0.98)

# --- Figure layout ---
using CairoMakie
fig = Figure(size=(1000, 760))

n = Observable(1)
title_str = @lift "t = " * prettytime(times[$n])
Label(fig[0, :], title_str, fontsize=18)

# Row 1: vertical cross-section (xz) at mid-domain y
ax_ωy = Axis(fig[1, 1]; title="Vorticity  ωy = ∂ᵤu − ∂ₓw  (xz slice)", xlabel="x", ylabel="z", aspect=DataAspect())
ax_w  = Axis(fig[1, 3]; title="Vertical velocity  w  (xz slice)", xlabel="x", ylabel="z", aspect=DataAspect())

ωy_plt = @lift view(ωy_ts[$n], :, Ny÷2, :)
w_plt  = @lift view(w_ts[$n], :, Ny÷2, :)

hm_ωy = heatmap!(ax_ωy, ωy_plt; colormap=:vik, colorrange=(-ωy_lim, ωy_lim))
hm_w  = heatmap!(ax_w,  w_plt;  colormap=:balance, colorrange=(-w_lim,  w_lim))

Colorbar(fig[1, 2], hm_ωy; label="ωy (s⁻¹)", vertical=true)
Colorbar(fig[1, 4], hm_w;  label="w (m s⁻¹)", vertical=true)

# Row 2: surface (xy) vorticity
ax_ωz = Axis(fig[2, 1]; title="Surface vorticity  ωz = ∂ₓv − ∂ᵧu", xlabel="x", ylabel="y", aspect=DataAspect())

ωz_plt = @lift view(ωz_ts[$n], :, :, Nz)

hm_ωz = heatmap!(ax_ωz, ωz_plt; colormap=:vik, colorrange=(-ωz_lim, ωz_lim))

Colorbar(fig[2, 2], hm_ωz; label="ωz (s⁻¹)", vertical=true)

# --- Record animation ---
record(fig, "hill_flow.mp4", 1:Nt; framerate=20) do nn
    n[] = nn
end

@info "Animation saved to hill_flow.mp4"
