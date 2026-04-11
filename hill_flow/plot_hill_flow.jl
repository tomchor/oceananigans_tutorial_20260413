using Oceananigans
import NCDatasets
using Statistics: quantile
using CairoMakie

# =============================================================================
# Animate hill_flow.jl output.
# Run hill_flow.jl first to produce hill_flow.nc.
# =============================================================================

ω_ts = FieldTimeSeries("hill_flow.nc", "ω")
w_ts = FieldTimeSeries("hill_flow.nc", "w")

times = ω_ts.times
Nt    = length(times)

ω_lim = max(quantile(abs.(vec(interior(ω_ts, :, 1, :, :))), 0.98), eps())
w_lim = max(quantile(abs.(vec(interior(w_ts, :, 1, :, :))), 0.98), eps())

# --- Figure layout ---
fig = Figure(size=(1000, 480))

n = Observable(1)
title_str = @lift "t = " * prettytime(times[$n])
Label(fig[0, 1:4], title_str, fontsize=18)

ax_ω = Axis(fig[1, 1]; title="Vorticity  ω = ∂ᵤu − ∂ₓw", xlabel="x", ylabel="z", aspect=DataAspect())
ax_w = Axis(fig[1, 3]; title="Vertical velocity  w",        xlabel="x", ylabel="z", aspect=DataAspect())

ω_plt = @lift view(ω_ts[$n], :, 1, :)
w_plt = @lift view(w_ts[$n], :, 1, :)

hm_ω = heatmap!(ax_ω, ω_plt; colormap=:vik,     colorrange=(-ω_lim, ω_lim))
hm_w = heatmap!(ax_w, w_plt; colormap=:balance,  colorrange=(-w_lim, w_lim))

Colorbar(fig[1, 2], hm_ω; label="ω (s⁻¹)",   vertical=true)
Colorbar(fig[1, 4], hm_w; label="w (m s⁻¹)", vertical=true)

# --- Record animation ---
record(fig, "hill_flow.mp4", 1:Nt; framerate=20) do nn
    n[] = nn
end

@info "Animation saved to hill_flow.mp4"
