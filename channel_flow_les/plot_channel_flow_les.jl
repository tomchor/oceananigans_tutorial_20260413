using Oceananigans
using CairoMakie
using Statistics

# =============================================================================
# Animate surface vorticity and tracer from channel_flow_les.jl output.
# Run channel_flow_les.jl first to produce channel_flow_surface.nc.
# =============================================================================

ζ_ts = FieldTimeSeries("channel_flow_surface.nc", "ζ")
c_ts = FieldTimeSeries("channel_flow_surface.nc", "c")

x, y, _ = nodes(ζ_ts)
times    = ζ_ts.times
Nt       = length(times)

# Symmetric color range for vorticity (computed once from all frames)
ζ_all = interior(ζ_ts, :, :, 1, :)
ζ_lim = quantile(abs.(vec(ζ_all)), 0.98)

# --- Figure layout ---
fig = Figure(size=(1000, 480))

ax_ζ = Axis(fig[1, 1]; title="Vertical vorticity  ζ = ∂ₓv − ∂ᵧu",
            xlabel="x", ylabel="y", aspect=DataAspect())
ax_c = Axis(fig[1, 3]; title="Tracer concentration  c",
            xlabel="x", ylabel="y", aspect=DataAspect())

n = Observable(1)

ζ_plt    = @lift interior(ζ_ts[$n], :, :, 1)
c_plt    = @lift interior(c_ts[$n], :, :, 1)
title_str = @lift "t = " * prettytime(times[$n])
Label(fig[0, :], title_str, fontsize=18)

hm_ζ = heatmap!(ax_ζ, x, y, ζ_plt; colormap=:vik,    colorrange=(-ζ_lim, ζ_lim))
hm_c = heatmap!(ax_c, x, y, c_plt; colormap=:thermal, colorrange=(0, 1))

Colorbar(fig[1, 2], hm_ζ; label="ζ (s⁻¹)", vertical=true)
Colorbar(fig[1, 4], hm_c; label="c",         vertical=true)

# --- Record animation ---
record(fig, "channel_flow_les.mp4", 1:Nt; framerate=20) do nn
    n[] = nn
end

@info "Animation saved to channel_flow_les.mp4"
