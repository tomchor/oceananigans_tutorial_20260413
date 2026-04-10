using Oceananigans
using CairoMakie
using Statistics: quantile

# =============================================================================
# Animate channel_flow_les.jl output.
# Run channel_flow_les.jl first to produce the .nc files.
# =============================================================================

# --- Surface fields ---
ζ_ts = FieldTimeSeries("channel_flow_surface.nc", "ζ")
c_ts = FieldTimeSeries("channel_flow_surface.nc", "c")

x_xy, y_xy, _ = nodes(ζ_ts)
times = ζ_ts.times
Nt    = length(times)

ζ_lim = quantile(abs.(vec(interior(ζ_ts, :, :, 1, :))), 0.98)

# --- xz slice ---
w_ts    = FieldTimeSeries("channel_flow_xz.nc", "w")
c_xz_ts = FieldTimeSeries("channel_flow_xz.nc", "c")

x_xz, _, z_xz = nodes(w_ts)

w_lim = quantile(abs.(vec(interior(w_ts, :, 1, :, :))), 0.98)

# --- Figure layout ---
fig = Figure(size=(1000, 760))

n = Observable(1)
title_str = @lift "t = " * prettytime(times[$n])
Label(fig[0, :], title_str, fontsize=18)

# Row 1: vertical cross-section (xz) at mid-domain y
ax_w    = Axis(fig[1, 1]; title="Vertical velocity  w  (xz slice)", xlabel="x", ylabel="z", aspect=DataAspect())
ax_cxz  = Axis(fig[1, 3]; title="Tracer  c  (xz slice)", xlabel="x", ylabel="z", aspect=DataAspect())

w_plt    = @lift w_ts[$n]
c_xz_plt = @lift c_xz_ts[$n]

hm_w   = heatmap!(ax_w,   w_plt;    colormap=:balance, colorrange=(-w_lim, w_lim))
hm_cxz = heatmap!(ax_cxz, c_xz_plt; colormap=:thermal, colorrange=(0, 1))

Colorbar(fig[1, 2], hm_w;   label="w (m s⁻¹)", vertical=true)
Colorbar(fig[1, 4], hm_cxz; label="c",          vertical=true)

# Row 2: surface (xy) fields
ax_ζ = Axis(fig[2, 1]; title="Surface vorticity  ζ = ∂ₓv − ∂ᵧu", xlabel="x", ylabel="y", aspect=DataAspect())
ax_c  = Axis(fig[2, 3]; title="Surface tracer  c", xlabel="x", ylabel="y", aspect=DataAspect())

ζ_plt = @lift ζ_ts[$n]
c_plt = @lift c_ts[$n]

hm_ζ = heatmap!(ax_ζ, ζ_plt; colormap=:vik,    colorrange=(-ζ_lim, ζ_lim))
hm_c  = heatmap!(ax_c, c_plt; colormap=:thermal, colorrange=(0, 1))

Colorbar(fig[2, 2], hm_ζ; label="ζ (s⁻¹)", vertical=true)
Colorbar(fig[2, 4], hm_c;  label="c",        vertical=true)

# --- Record animation ---
record(fig, "channel_flow_les.mp4", 1:Nt; framerate=20) do nn
    n[] = nn
end

@info "Animation saved to channel_flow_les.mp4"
