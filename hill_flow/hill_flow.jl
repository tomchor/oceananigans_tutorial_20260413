using Oceananigans
using Oceananigans.Solvers: ConjugateGradientPoissonSolver
using NCDatasets
using Printf

# =============================================================================
# Flow past a Gaussian hill
#
# A 3-D nonhydrostatic simulation with open east/west boundaries.
# A background flow U∞ enters from the west, goes over a hill, and exits
# to the east.  The immersed boundary method handles the topography.
#
# Stratification (N² > 0) can be added to show internal wave generation.
# =============================================================================

# --- Physical parameters ---
Lx = 20.0
Ly = Lx/2
H  = 2.0
U∞ = 1.0
z₀ = 1e-4       # roughness length (m)

# Hill geometry (axisymmetric Gaussian, centered in the domain)
x₀ = 0.0        # x center position
h₀ = 0.6H       # peak height above the bottom
σ  = Lx / 10    # horizontal half-width

hill(x, y, p) = p.h₀ * exp(-((x - p.x₀)^2 + y^2) / p.σ^2) - p.H   # returns z_bottom(x, y)

# --- Grid ---
Nx, Ny, Nz = 128, 128, 32

underlying_grid = RectilinearGrid(size     = (Nx, Ny, Nz),
                                  x        = (-Lx/2, Lx/2),
                                  y        = (-Ly/2, Ly/2),
                                  z        = (-H, 0),
                                  topology = (Bounded, Periodic, Bounded),
                                  halo     = (6, 6, 6))

# Drag coefficient from law of the wall (https://doi.org/10.1029/2005WR004685)
const κᵛᵏ = 0.4    # von Kármán constant
z₁ = minimum_zspacing(underlying_grid, Center(), Center(), Center()) / 2
Cd = (κᵛᵏ / log(z₁ / z₀))^2
@info "z₁ = $z₁,  Cd = $Cd"

hill_params = (; x₀, h₀, σ, H)
grid = ImmersedBoundaryGrid(underlying_grid, GridFittedBottom((x, y) -> hill(x, y, hill_params)))

# --- Boundary conditions ---
drag  = BulkDrag(coefficient=Cd)
u_bcs = FieldBoundaryConditions(west   = OpenBoundaryCondition(U∞),
                                east   = OpenBoundaryCondition(U∞, scheme = PerturbationAdvection()),
                                bottom = drag)
v_bcs = FieldBoundaryConditions(bottom = drag)

# --- Model ---
model = NonhydrostaticModel(grid;
                            #pressure_solver= ConjugateGradientPoissonSolver(grid; maxiter = 10),
                            boundary_conditions = (u=u_bcs, v=v_bcs),
                            advection           = WENO(order=5),
                            timestepper         = :RungeKutta3)

set!(model, u=U∞)

# --- Simulation ---
Δt₀ = 0.1 * minimum_xspacing(grid) / U∞
simulation = Simulation(model; Δt=Δt₀, stop_time=30)
conjure_time_step_wizard!(simulation, cfl=0.8, IterationInterval(5))

wall_clock = Ref(time_ns())
function progress(sim)
    u = sim.model.velocities.u
    elapsed = prettytime(1e-9 * (time_ns() - wall_clock[]))
    @info @sprintf("t = %s, Δt = %s, max|u| = %.3f, elapsed wall time = %s",
                   prettytime(time(sim)), prettytime(sim.Δt), maximum(abs, u), elapsed)
    wall_clock[] = time_ns()
end
add_callback!(simulation, progress, IterationInterval(100))

# --- Output ---
u, v, w = model.velocities
ωy = Field(∂z(u) - ∂x(w))    # y-component of vorticity (in the xz plane)
ωz = Field(∂x(v) - ∂y(u))    # z-component of vorticity (in the xy plane)

simulation.output_writers[:fields] = NetCDFWriter(model,
    (; u, v, w, ωy, ωz),
    schedule           = TimeInterval(0.5),
    filename           = "hill_flow.nc",
    overwrite_existing = true)

run!(simulation)
