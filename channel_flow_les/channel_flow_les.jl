using Oceananigans
using NCDatasets
using Printf

# =============================================================================
# Forced channel flow with quadratic bottom drag
#
# A doubly-periodic LES channel driven by a constant pressure-gradient body
# force and damped by a quadratic bottom drag.  A passive tracer blob is
# released at the surface so you can watch it get stirred by the turbulence.
#
# Increase Nx/Ny/Nz for better-resolved turbulence.
# =============================================================================

# --- Physical parameters ---
const Lx = 2π          # m, domain length (zonal)
const Ly = 2π          # m, domain length (meridional)
const H  = 1.0         # m, depth
const U₀ = 1.0         # m/s, target mean velocity
const Cd = 1e-2        # quadratic drag coefficient

# Equilibrium forcing: balances drag at mean flow U₀  (F₀ = Cd U₀² / H)
const F₀ = Cd * U₀^2 / H

# --- Grid ---
Nx, Ny, Nz = 64, 64, 32

grid = RectilinearGrid(size  = (Nx, Ny, Nz),
                       x     = (0, Lx),
                       y     = (0, Ly),
                       z     = (-H, 0),
                       topology = (Periodic, Periodic, Bounded))

# --- Quadratic bottom drag (applied as a bottom flux) ---
@inline drag_u(x, y, t, u, v) = -Cd * u * sqrt(u^2 + v^2)
@inline drag_v(x, y, t, u, v) = -Cd * v * sqrt(u^2 + v^2)

u_bcs = FieldBoundaryConditions(bottom = FluxBoundaryCondition(drag_u, field_dependencies=(:u, :v)))
v_bcs = FieldBoundaryConditions(bottom = FluxBoundaryCondition(drag_v, field_dependencies=(:u, :v)))

# --- Pressure-gradient body force (drives flow in x) ---
@inline pressure_forcing(x, y, z, t) = F₀

# --- Model ---
model = NonhydrostaticModel(; grid,
                              closure              = SmagorinskyLilly(),
                              advection            = WENO(order=5),
                              timestepper          = :RungeKutta3,
                              tracers              = :c,
                              boundary_conditions  = (u=u_bcs, v=v_bcs),
                              forcing              = (u=Forcing(pressure_forcing),))

# --- Initial conditions ---
uᵢ(x, y, z) = U₀ * (1 + 0.05 * randn())
cᵢ(x, y, z) = exp(-((x - Lx/4)^2 + (y - Ly/4)^2) / (Lx/8)^2)  # Gaussian tracer blob
set!(model, u=uᵢ, c=cᵢ)

# --- Simulation ---
simulation = Simulation(model; Δt=1e-2, stop_time=50)
conjure_time_step_wizard!(simulation, cfl=0.8, IterationInterval(5))

wall_clock = Ref(time_ns())
function progress(sim)
    u = sim.model.velocities.u
    elapsed = prettytime(1e-9 * (time_ns() - wall_clock[]))
    @info @sprintf("t = %s, Δt = %s, max|u| = %.3f, wall time = %s",
                   prettytime(time(sim)), prettytime(sim.Δt), maximum(abs, u), elapsed)
    wall_clock[] = time_ns()
end
add_callback!(simulation, progress, IterationInterval(100))

# --- Output: surface (xy) slice ---
u, v, w = model.velocities
ζ = Field(∂x(v) - ∂y(u))   # vertical vorticity at the surface

simulation.output_writers[:surface] = NetCDFWriter(model,
    merge(model.velocities, model.tracers, (; ζ)),
    schedule         = TimeInterval(1.0),
    filename         = "channel_flow_surface.nc",
    indices          = (:, :, Nz),       # save only the surface layer
    overwrite_existing = true)

run!(simulation)
