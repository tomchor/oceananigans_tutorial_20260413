using Oceananigans
import NCDatasets

# =============================================================================
# Animate hill_flow.jl output.
# Run hill_flow.jl first to produce hill_flow.nc.
# =============================================================================

ω_ts = FieldTimeSeries("hill_flow.nc", "ω")
w_ts = FieldTimeSeries("hill_flow.nc", "w")
u_ts = FieldTimeSeries("hill_flow.nc", "u")

times = ω_ts.times
Nt    = length(times)

using Statistics: quantile
ω_lim = max(quantile(abs.(vec(interior(ω_ts, :, 1, :, :))), 0.98), eps())
w_lim = max(quantile(abs.(vec(interior(w_ts, :, 1, :, :))), 0.98), eps())
u_lim = max(quantile(abs.(vec(interior(u_ts, :, 1, :, :))), 0.98), eps())

# --- Figure layout ---
using GLMakie
fig = Figure(size=(900, 500))
n = Observable(1)

using Printf
title_str = @lift @sprintf("t = %.1f", times[$n])
Label(fig[0, 1:2], title_str, fontsize=18)

ax_ω = Axis(fig[1, 1]; title="Vorticity  ω = ∂ᵤu − ∂ₓw", xlabel="x", ylabel="z", aspect=DataAspect())
ax_w = Axis(fig[2, 1]; title="Vertical velocity  w",        xlabel="x", ylabel="z", aspect=DataAspect())
ax_u = Axis(fig[3, 1]; title="Horizontal velocity  u",      xlabel="x", ylabel="z", aspect=DataAspect())

ω_plt = @lift view(ω_ts[$n], :, 1, :)
w_plt = @lift view(w_ts[$n], :, 1, :)
u_plt = @lift view(u_ts[$n], :, 1, :)

hm_ω = heatmap!(ax_ω, ω_plt; colormap=:vik,     colorrange=(-ω_lim, ω_lim))
hm_w = heatmap!(ax_w, w_plt; colormap=:balance,  colorrange=(-w_lim, w_lim))
hm_u = heatmap!(ax_u, u_plt; colormap=:thermal,  colorrange=(0, u_lim))

Colorbar(fig[1, 2], hm_ω; label="ω", vertical=true, height=Relative(0.3))
Colorbar(fig[2, 2], hm_w; label="w", vertical=true, height=Relative(0.3))
Colorbar(fig[3, 2], hm_u; label="u", vertical=true, height=Relative(0.3))

# --- Record animation ---
record(fig, "hill_flow.mp4", 1:Nt; framerate=20) do nn
    n[] = nn
end

@info "Animation saved to hill_flow.mp4"
