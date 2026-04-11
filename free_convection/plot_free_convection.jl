using CairoMakie
using Printf
using Oceananigans
using Statistics: quantile

# =============================================================================
# Animate free_convection.jl output.
# Run free_convection.jl first to produce free_convection.nc and
# free_convection_xz.nc.
# =============================================================================

# --- Load timeseries ---
w_xy = FieldTimeSeries("free_convection.nc",    "w")
T_xy = FieldTimeSeries("free_convection.nc",    "T")
w_xz = FieldTimeSeries("free_convection_xz.nc", "w")
T_xz = FieldTimeSeries("free_convection_xz.nc", "T")

x, y, _ = nodes(w_xy)
x2, _, z = nodes(w_xz)
times = w_xy.times
Nt    = length(times)

# Compute shared colorranges across both slices
w_all = vcat(vec(interior(w_xy, :, :, 1, :)), vec(interior(w_xz, :, 1, :, :)))
T_all = vcat(vec(interior(T_xy, :, :, 1, :)), vec(interior(T_xz, :, 1, :, :)))

w_lim = quantile(abs.(w_all), 0.98)
T_min = quantile(T_all, 0.02)
T_max = quantile(T_all, 0.98)

# --- Figure layout ---
fig = Figure(size=(1000, 860))

n = Observable(1)
title_str = @lift "t = " * prettytime(times[$n])
Label(fig[0, 1:4], title_str, fontsize=18)

kwargs_xy = (xlabel="x (m)", ylabel="y (m)", aspect=DataAspect())
kwargs_xz = (xlabel="x (m)", ylabel="z (m)", aspect=DataAspect())

ax_w_xy = Axis(fig[1, 1]; title="w  —  horizontal slice (z ≈ H/2)", kwargs_xy...)
ax_T_xy = Axis(fig[1, 3]; title="θ  —  horizontal slice (z ≈ H/2)", kwargs_xy...)
ax_w_xz = Axis(fig[2, 1]; title="w  —  vertical cross-section",     kwargs_xz...)
ax_T_xz = Axis(fig[2, 3]; title="θ  —  vertical cross-section",     kwargs_xz...)

w_xy_plt = @lift interior(w_xy[$n], :, :, 1)
T_xy_plt = @lift interior(T_xy[$n], :, :, 1)
w_xz_plt = @lift interior(w_xz[$n], :, 1, :)
T_xz_plt = @lift interior(T_xz[$n], :, 1, :)

hm_w_xy = heatmap!(ax_w_xy, x,  y, w_xy_plt; colormap=:vik,     colorrange=(-w_lim, w_lim))
hm_T_xy = heatmap!(ax_T_xy, x,  y, T_xy_plt; colormap=:thermal,  colorrange=(T_min, T_max))
hm_w_xz = heatmap!(ax_w_xz, x2, z, w_xz_plt; colormap=:vik,     colorrange=(-w_lim, w_lim))
hm_T_xz = heatmap!(ax_T_xz, x2, z, T_xz_plt; colormap=:thermal,  colorrange=(T_min, T_max))

Colorbar(fig[1, 2], hm_w_xy; label="w (m/s)", vertical=true)
Colorbar(fig[1, 4], hm_T_xy; label="θ (K)",   vertical=true)
Colorbar(fig[2, 2], hm_w_xz; label="w (m/s)", vertical=true)
Colorbar(fig[2, 4], hm_T_xz; label="θ (K)",   vertical=true)

# --- Record animation ---
record(fig, "free_convection.mp4", 1:Nt; framerate=15) do nn
    n[] = nn
end

@info "Animation saved to free_convection.mp4"
