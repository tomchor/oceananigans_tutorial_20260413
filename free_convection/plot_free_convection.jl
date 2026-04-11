using CairoMakie
using Printf
using Oceananigans
using Statistics: quantile

# =============================================================================
# Animate free_convection.jl output.
# Run free_convection.jl first to produce free_convection.nc.
# =============================================================================

w_ts = FieldTimeSeries("free_convection.nc", "w")
T_ts = FieldTimeSeries("free_convection.nc", "T")

x, y, _ = nodes(w_ts)
times    = w_ts.times
Nt       = length(times)

w_lim = quantile(abs.(vec(interior(w_ts, :, :, 1, :))), 0.98)
T_min = quantile(vec(interior(T_ts, :, :, 1, :)), 0.02)
T_max = quantile(vec(interior(T_ts, :, :, 1, :)), 0.98)

# --- Figure layout ---
fig = Figure(size=(1000, 480))

n = Observable(1)
title_str = @lift "t = " * prettytime(times[$n])
Label(fig[0, :], title_str, fontsize=18)

ax_w = Axis(fig[1, 1]; title="Vertical velocity  w  (m/s)",
            xlabel="x (m)", ylabel="y (m)", aspect=DataAspect())
ax_T = Axis(fig[1, 3]; title="Potential temperature  θ  (K)",
            xlabel="x (m)", ylabel="y (m)", aspect=DataAspect())

w_plt = @lift interior(w_ts[$n], :, :, 1)
T_plt = @lift interior(T_ts[$n], :, :, 1)

hm_w = heatmap!(ax_w, x, y, w_plt; colormap=:vik,    colorrange=(-w_lim, w_lim))
hm_T = heatmap!(ax_T, x, y, T_plt; colormap=:thermal, colorrange=(T_min, T_max))

Colorbar(fig[1, 2], hm_w; label="w (m/s)", vertical=true)
Colorbar(fig[1, 4], hm_T; label="θ (K)",   vertical=true)

# --- Record animation ---
record(fig, "free_convection.mp4", 1:Nt; framerate=15) do nn
    n[] = nn
end

@info "Animation saved to free_convection.mp4"
