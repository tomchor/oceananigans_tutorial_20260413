using CairoMakie
using NCDatasets

# =============================================================================
# Animate surface vorticity and tracer from channel_flow_les.jl output.
# Run channel_flow_les.jl first to produce channel_flow_surface.nc.
# =============================================================================

ds = NCDataset("channel_flow_surface.nc")

x = ds["xC"][:]   # cell-center x positions
y = ds["yC"][:]   # cell-center y positions
t = ds["time"][:]

# Variables are saved as (x, y, z=1, time); squeeze the singleton z dimension.
ζ = dropdims(ds["ζ"][:, :, 1, :], dims=3)   # (Nx, Ny, Nt)  vertical vorticity
c = dropdims(ds["c"][:, :, 1, :], dims=3)   # (Nx, Ny, Nt)  passive tracer

close(ds)

Nt = length(t)

# Symmetric color range for vorticity (computed once from all frames)
ζ_lim = quantile(abs.(vec(ζ)), 0.98)

# --- Figure layout ---
fig = Figure(size=(1000, 480))

ax_ζ = Axis(fig[1, 1]; title="Vertical vorticity  ζ = ∂ₓv − ∂ᵧu",
            xlabel="x", ylabel="y", aspect=DataAspect())
ax_c = Axis(fig[1, 3]; title="Tracer concentration  c",
            xlabel="x", ylabel="y", aspect=DataAspect())

n = Observable(1)

ζ_plt = @lift ζ[:, :, $n]
c_plt = @lift c[:, :, $n]
title_str = @lift @sprintf("t = %.1f s", t[$n])
Label(fig[0, :], title_str, fontsize=18)

hm_ζ = heatmap!(ax_ζ, x, y, ζ_plt; colormap=:vik,   colorrange=(-ζ_lim, ζ_lim))
hm_c = heatmap!(ax_c, x, y, c_plt; colormap=:thermal, colorrange=(0, 1))

Colorbar(fig[1, 2], hm_ζ; label="ζ (s⁻¹)", vertical=true)
Colorbar(fig[1, 4], hm_c; label="c",         vertical=true)

# --- Record animation ---
record(fig, "channel_flow_les.mp4", 1:Nt; framerate=20) do nn
    n[] = nn
end

@info "Animation saved to channel_flow_les.mp4"
