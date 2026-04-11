using CairoMakie
using Printf
using Oceananigans
using Statistics: quantile

# =============================================================================
# Animate kelvin_helmholtz.jl output.
# Run kelvin_helmholtz.jl first to produce kelvin_helmholtz.nc.
# =============================================================================

plot_filepath = "kelvin_helmholtz.nc"

#+++ Load timeseries
@info "Loading timeseries..."
ω_timeseries = FieldTimeSeries(plot_filepath, "ω")
b_timeseries = FieldTimeSeries(plot_filepath, "b")
S_timeseries = FieldTimeSeries(plot_filepath, "S")

times = ω_timeseries.times
S_lim = quantile(vec(interior(S_timeseries, :, 1, :, :)), 0.98)
#---

#+++ Build figure
n = Observable(1)

ωₙ = @lift view(ω_timeseries[$n], :, 1, :)
bₙ = @lift view(b_timeseries[$n], :, 1, :)
Sₙ = @lift view(S_timeseries[$n], :, 1, :)

fig = Figure(size=(1200, 500))

title = @lift @sprintf("Kelvin-Helmholtz Instability\nt = %.1f", times[$n])
fig[1, 1:6] = Label(fig, title, fontsize=20, tellwidth=false, justification=:center)

kwargs = (xlabel="x", ylabel="z", aspect=1)

ax_ω = Axis(fig[2, 1]; title="Vorticity",      kwargs...)
ax_b = Axis(fig[2, 3]; title="Buoyancy",        kwargs...)
ax_S = Axis(fig[2, 5]; title="Strain rate (S)", kwargs...)

hm_ω = heatmap!(ax_ω, ωₙ; colormap=:balance, colorrange=(-1, 1))
Colorbar(fig[2, 2], hm_ω)

hm_b = heatmap!(ax_b, bₙ; colormap=:balance, colorrange=(-0.08, 0.08))
Colorbar(fig[2, 4], hm_b)

hm_S = heatmap!(ax_S, Sₙ; colormap=:thermal, colorrange=(0, S_lim))
Colorbar(fig[2, 6], hm_S)
#---

#+++ Record animation
animation_filename = "kelvin_helmholtz.mp4"
record(fig, animation_filename, 1:length(times); framerate=12) do i
    @info "Plotting frame $i of $(length(times))..."
    n[] = i
end
@info "Animation saved as $(animation_filename)"
