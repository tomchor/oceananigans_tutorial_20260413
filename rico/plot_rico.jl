using Oceananigans
using GLMakie
using Statistics: quantile
using Printf

# =============================================================================
# Visualize rico.jl output.
# Run rico.jl first to produce rico.jld2.
# Produces:
#   rico.mp4 – animation of xz cross-sections
# =============================================================================

plot_filepath = "rico.jld2"

# --- Load timeseries ---
@info "Loading timeseries..."
w_ts   = FieldTimeSeries(plot_filepath, "w")
θ_ts   = FieldTimeSeries(plot_filepath, "θ")
qᶜˡ_ts = FieldTimeSeries(plot_filepath, "qᶜˡ")
qʳ_ts  = FieldTimeSeries(plot_filepath, "qʳ")

times = θ_ts.times
Nt    = length(times)

# --- Colormap limits ---
w_lim   = quantile(abs.(vec(interior(w_ts))), 0.998)
qᶜˡ_lim = quantile(vec(interior(qᶜˡ_ts)), 0.999)
qʳ_lim  = quantile(vec(interior(qʳ_ts)), 0.99)
θ_min, θ_max = extrema(interior(θ_ts))

# =============================================================================
# Part 1: animation of xz fields
# =============================================================================
@info "Building animation figure..."

n = Observable(1)
title_str = @lift @sprintf("RICO  —  t = %s", prettytime(times[$n]))

w_n   = @lift w_ts[$n]
θ_n   = @lift θ_ts[$n]
qᶜˡ_n = @lift qᶜˡ_ts[$n]
qʳ_n  = @lift qʳ_ts[$n]

fig = Figure(size=(1200, 700))
Label(fig[0, :], title_str, fontsize=18, tellwidth=false)

kw = (xlabel="x (m)", ylabel="z (m)")

ax_w   = Axis(fig[1, 1]; title="Vertical velocity  w",   kw...)
ax_θ   = Axis(fig[1, 3]; title="Potential temperature θ", kw...)
ax_qᶜˡ = Axis(fig[2, 1]; title="Cloud liquid  qᶜˡ",     kw...)
ax_qʳ  = Axis(fig[2, 3]; title="Rain  qʳ",              kw...)

hm_w   = heatmap!(ax_w,   w_n;   colormap=:balance, colorrange=(-w_lim, w_lim))
hm_θ   = heatmap!(ax_θ,   θ_n;   colormap=:thermal, colorrange=(θ_min, θ_max))
hm_qᶜˡ = heatmap!(ax_qᶜˡ, qᶜˡ_n; colormap=:dense, colorrange=(0, qᶜˡ_lim))
hm_qʳ  = heatmap!(ax_qʳ,  qʳ_n;  colormap=:amp, colorrange=(0, qʳ_lim))

Colorbar(fig[1, 2], hm_w;   label="w (m s⁻¹)")
Colorbar(fig[1, 4], hm_θ;   label="θ (K)")
Colorbar(fig[2, 2], hm_qᶜˡ; label="qᶜˡ (kg kg⁻¹)")
Colorbar(fig[2, 4], hm_qʳ;  label="qʳ (kg kg⁻¹)")

record(fig, "rico.mp4", 1:Nt; framerate=12) do i
    @info @sprintf("Animating frame %d / %d  (t = %s)", i, Nt, prettytime(times[i]))
    n[] = i
end
@info "Animation saved as rico.mp4"
