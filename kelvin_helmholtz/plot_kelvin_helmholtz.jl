using Oceananigans
using CairoMakie
using Statistics: quantile

# =============================================================================
# Animate kelvin_helmholtz.jl output.
# Run kelvin_helmholtz.jl first to produce kelvin_helmholtz.nc.
# =============================================================================

ω_ts = FieldTimeSeries("kelvin_helmholtz.nc", "ω")
b_ts = FieldTimeSeries("kelvin_helmholtz.nc", "b")
S_ts = FieldTimeSeries("kelvin_helmholtz.nc", "S")

x, _, z = nodes(ω_ts)
times = ω_ts.times
Nt    = length(times)

ω_lim = quantile(abs.(vec(interior(ω_ts, :, 1, :, :))), 0.98)
b_lim = maximum(abs.(vec(interior(b_ts, :, 1, :, :))))
S_lim = quantile(vec(interior(S_ts, :, 1, :, :)), 0.98)

# --- Figure layout ---
fig = Figure(size=(1400, 480))

n = Observable(1)
title_str = @lift "t = " * prettytime(times[$n])
Label(fig[0, :], title_str, fontsize=18)

ax_ω = Axis(fig[1, 1]; title="Vorticity  ω = ∂ᵤu − ∂ₓw",
            xlabel="x", ylabel="z", aspect=DataAspect())
ax_b = Axis(fig[1, 3]; title="Buoyancy  b",
            xlabel="x", ylabel="z", aspect=DataAspect())
ax_S = Axis(fig[1, 5]; title="Strain rate modulus  S",
            xlabel="x", ylabel="z", aspect=DataAspect())

ω_plt = @lift interior(ω_ts[$n], :, 1, :)
b_plt = @lift interior(b_ts[$n], :, 1, :)
S_plt = @lift interior(S_ts[$n], :, 1, :)

hm_ω = heatmap!(ax_ω, x, z, ω_plt; colormap=:vik,     colorrange=(-ω_lim, ω_lim))
hm_b = heatmap!(ax_b, x, z, b_plt; colormap=:balance,  colorrange=(-b_lim, b_lim))
hm_S = heatmap!(ax_S, x, z, S_plt; colormap=:thermal,  colorrange=(0, S_lim))

Colorbar(fig[1, 2], hm_ω; label="ω (s⁻¹)",   vertical=true)
Colorbar(fig[1, 4], hm_b; label="b (m s⁻²)", vertical=true)
Colorbar(fig[1, 6], hm_S; label="S (s⁻¹)",   vertical=true)

# --- Record animation ---
record(fig, "kelvin_helmholtz.mp4", 1:Nt; framerate=20) do nn
    n[] = nn
end

@info "Animation saved to kelvin_helmholtz.mp4"
