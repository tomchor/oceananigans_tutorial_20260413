using Oceananigans
using GLMakie
using Statistics: mean, quantile
using Printf

# =============================================================================
# Visualize rico.jl output.
# Run rico.jl first to produce rico.jld2.
# Produces:
#   rico.mp4         – animation of xz cross-sections
#   rico_profiles.png – time-evolving x-averaged profiles
# =============================================================================

plot_filepath = "rico.jld2"

# --- Load timeseries ---
@info "Loading timeseries..."
u_ts   = FieldTimeSeries(plot_filepath, "u")
w_ts   = FieldTimeSeries(plot_filepath, "w")
θ_ts   = FieldTimeSeries(plot_filepath, "θ")
qᶜˡ_ts = FieldTimeSeries(plot_filepath, "qᶜˡ")
qʳ_ts  = FieldTimeSeries(plot_filepath, "qʳ")

times = θ_ts.times
Nt    = length(times)

# --- Colormap limits ---
w_lim   = max(quantile(abs.(vec(interior(w_ts))),  0.998), eps())
qᶜˡ_lim = max(quantile(vec(interior(qᶜˡ_ts)),     0.99),  eps())
qʳ_lim  = max(quantile(vec(interior(qʳ_ts)),      0.99),  eps())
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
hm_qᶜˡ = heatmap!(ax_qᶜˡ, qᶜˡ_n; colormap=:dense,  colorrange=(0, qᶜˡ_lim))
hm_qʳ  = heatmap!(ax_qʳ,  qʳ_n;  colormap=:amp,    colorrange=(0, qʳ_lim))

Colorbar(fig[1, 2], hm_w;   label="w (m s⁻¹)")
Colorbar(fig[1, 4], hm_θ;   label="θ (K)")
Colorbar(fig[2, 2], hm_qᶜˡ; label="qᶜˡ (kg kg⁻¹)")
Colorbar(fig[2, 4], hm_qʳ;  label="qʳ (kg kg⁻¹)")

record(fig, "rico.mp4", 1:Nt; framerate=12) do i
    @info @sprintf("Animating frame %d / %d  (t = %s)", i, Nt, prettytime(times[i]))
    n[] = i
end
@info "Animation saved as rico.mp4"

# =============================================================================
# Part 2: time-evolving x-averaged profiles
# =============================================================================
@info "Building profiles figure..."

z = Oceananigans.Grids.znodes(θ_ts.grid, Center())

# Select ~8 evenly-spaced snapshots
step    = max(1, Nt ÷ 8)
indices = 1:step:Nt
colors  = cgrad(:viridis, length(indices); categorical=true)

fig2 = Figure(size=(1100, 480))
Label(fig2[0, :], "RICO: x-averaged profiles", fontsize=18, tellwidth=false)

ax_θ   = Axis(fig2[1, 1]; xlabel="θ (K)",           ylabel="z (m)")
ax_u   = Axis(fig2[1, 2]; xlabel="u (m s⁻¹)",       ylabel="z (m)")
ax_qᶜˡ = Axis(fig2[1, 3]; xlabel="qᶜˡ (kg kg⁻¹)",   ylabel="z (m)")
ax_qʳ  = Axis(fig2[1, 4]; xlabel="qʳ (kg kg⁻¹)",    ylabel="z (m)")

for (ci, ni) in enumerate(indices)
    label = prettytime(times[ni])
    c     = colors[ci]

    θ_prof   = dropdims(mean(interior(θ_ts[ni]),   dims=(1, 2)), dims=(1, 2))
    u_prof   = dropdims(mean(interior(u_ts[ni]),   dims=(1, 2)), dims=(1, 2))
    qᶜˡ_prof = dropdims(mean(interior(qᶜˡ_ts[ni]), dims=(1, 2)), dims=(1, 2))
    qʳ_prof  = dropdims(mean(interior(qʳ_ts[ni]),  dims=(1, 2)), dims=(1, 2))

    lines!(ax_θ,   θ_prof,   z; color=c, label=label)
    lines!(ax_u,   u_prof,   z; color=c)
    lines!(ax_qᶜˡ, qᶜˡ_prof, z; color=c)
    lines!(ax_qʳ,  qʳ_prof,  z; color=c)
end

for ax in (ax_θ, ax_u, ax_qᶜˡ, ax_qʳ)
    ylims!(ax, 0, 4000)
end

axislegend(ax_θ, position=:rb, labelsize=11)

save("rico_profiles.png", fig2)
@info "Profiles saved as rico_profiles.png"
