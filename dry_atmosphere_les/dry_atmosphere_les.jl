using Oceananigans
using Oceananigans.Units
using NCDatasets
using Printf

# =============================================================================
# Atmosphere LES with quadratic surface drag
#
# A doubly-periodic LES of the atmospheric boundary layer driven by a
# constant pressure-gradient body force and damped by a quadratic surface drag.
# A passive tracer blob is released near the surface so you can watch it get
# stirred by the turbulence.
#
# Increase Nx/Ny/Nz for better-resolved turbulence.
# =============================================================================

# --- Physical parameters ---
Lx = 2kilometers   # domain length (zonal)
Ly = 1kilometers   # domain length (meridional)
H  = 1kilometer    # boundary layer height
U₀ = 8.0          # m/s, target mean wind speed
z₀ = 0.1meters    # aerodynamic roughness length (moderate terrain / farmland)

# --- Grid ---
Nx, Ny, Nz = 64, 32, 32

grid = RectilinearGrid(size     = (Nx, Ny, Nz),
                       x        = (0, Lx),
                       y        = (0, Ly),
                       z        = (0, H),
                       topology = (Periodic, Periodic, Bounded))

# Drag coefficient from law of the wall (https://doi.org/10.1029/2005WR004685)
const κᵛᵏ = 0.4    # von Kármán constant
z₁ = minimum_zspacing(grid, Center(), Center(), Center()) / 2
Cd = (κᵛᵏ / log(z₁ / z₀))^2
@info "z₁ = $(z₁) m,  Cd = $(round(Cd, sigdigits=3))"

F₀ = Cd * U₀^2 / H   # pressure-gradient force balancing drag at mean wind U₀

# --- Quadratic surface drag ---
drag  = BulkDrag(coefficient=Cd)
u_bcs = FieldBoundaryConditions(bottom = drag)
v_bcs = FieldBoundaryConditions(bottom = drag)

# --- Pressure-gradient body force (drives flow in x) ---
@inline pressure_forcing(x, y, z, t, p) = p.F₀

# --- Model ---
model = NonhydrostaticModel(grid;
                            closure             = SmagorinskyLilly(),
                            advection           = WENO(order=5),
                            timestepper         = :RungeKutta3,
                            tracers             = :c,
                            boundary_conditions = (u=u_bcs, v=v_bcs),
                            forcing             = (u=Forcing(pressure_forcing, parameters=(; F₀)),))

# --- Initial conditions ---
uᵢ(x, y, z) = U₀ * (1 + 0.05 * randn() * exp(-z / (H/4)))
cᵢ(x, y, z) = exp(-((x - Lx/4)^2 + (y - Ly/4)^2) / (Lx/8)^2)  # Gaussian tracer blob
set!(model, u=uᵢ, c=cᵢ)

# --- Simulation ---
simulation = Simulation(model; Δt=5seconds, stop_time=4hours)
conjure_time_step_wizard!(simulation, cfl=0.8, IterationInterval(5))

wall_clock = Ref(time_ns())
function progress(sim)
    u = sim.model.velocities.u
    elapsed = prettytime(1e-9 * (time_ns() - wall_clock[]))
    @info @sprintf("t = %s, Δt = %s, max|u| = %.2f m/s, wall time = %s",
                   prettytime(time(sim)), prettytime(sim.Δt), maximum(abs, u), elapsed)
    wall_clock[] = time_ns()
end
add_callback!(simulation, progress, IterationInterval(100))

# --- Output ---
u, v, w = model.velocities
ζ = Field(∂x(v) - ∂y(u))   # vertical vorticity

simulation.output_writers[:surface] = NetCDFWriter(model,
    merge(model.velocities, model.tracers, (; ζ)),
    schedule            = TimeInterval(10minutes),
    filename            = "dry_atmosphere_surface.nc",
    indices             = (:, :, 2),        # near-surface (xy) slice
    overwrite_existing  = true)

simulation.output_writers[:xz_slice] = NetCDFWriter(model,
    merge(model.velocities, model.tracers),
    schedule            = TimeInterval(10minutes),
    filename            = "dry_atmosphere_xz.nc",
    indices             = (:, Ny÷2, :),     # xz slice at mid-domain y
    overwrite_existing  = true)

run!(simulation)
