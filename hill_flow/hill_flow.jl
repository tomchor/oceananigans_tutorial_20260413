using Oceananigans
using NCDatasets
using Oceanostics.ProgressMessengers

# =============================================================================
# Flow past a Gaussian hill (2D, xz)
#
# A nonhydrostatic simulation with open east/west boundaries.
# A background flow U∞ enters from the west, passes over a hill, and exits
# to the east.  The immersed boundary method handles the topography.
# =============================================================================

# --- Physical parameters ---
Lx = 20.0
H  = 2.0
U∞ = 1.0
z₀ = 1e-4       # roughness length (m)

# Hill geometry (Gaussian ridge, centered in the domain)
x₀ = 0.0        # x center position
h₀ = 0.4H       # peak height above the bottom
σ  = Lx / 10    # horizontal half-width

hill(x) = h₀ * exp(-((x - x₀) / σ)^2) - H   # returns z_bottom(x)

# --- Grid ---
Nx, Nz = 128, 32

underlying_grid = RectilinearGrid(size     = (Nx, Nz),
                                  x        = (-Lx/3, 2Lx/3),
                                  z        = (-H, 0),
                                  topology = (Bounded, Flat, Bounded))

const κᵛᵏ = 0.4 # von Kármán constant
z₁ = minimum_zspacing(underlying_grid, Center(), Center(), Center()) / 2
Cd = (κᵛᵏ / log(z₁ / z₀))^2
@info "z₁ = $z₁,  Cd = $Cd"

hill_params = (; x₀, h₀, σ, H)
grid = ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(hill))

# --- Boundary conditions ---
drag = BulkDrag(coefficient=Cd)
u_bcs = FieldBoundaryConditions(west     = OpenBoundaryCondition(U∞), # Constant inflow
                                east     = OpenBoundaryCondition(U∞, scheme = PerturbationAdvection()), # Perturbations get advected out
                                bottom   = drag,
                                immersed = drag)

# --- Model ---
using Oceananigans.Solvers: ConjugateGradientPoissonSolver
model = NonhydrostaticModel(grid;
                            boundary_conditions = (u=u_bcs,),
                            advection           = WENO(order=5), # Implitict dissipation
                            timestepper         = :RungeKutta3,
                            pressure_solver     = ConjugateGradientPoissonSolver(grid; maxiter = 10)) # More accurate results with immersed boundary

set!(model, u=U∞)

# --- Simulation ---
Δt₀ = 0.1 * minimum_xspacing(grid) / U∞
simulation = Simulation(model; Δt=Δt₀, stop_time=50)
conjure_time_step_wizard!(simulation, cfl=0.5, IterationInterval(2))

progress = SimulationTime() + TimeStep() + Walltime() + AdvectiveCFLNumber() + MaxUVelocity()
add_callback!(simulation, progress, IterationInterval(100))

# --- Output ---
u, v, w = model.velocities
ω = Field(∂z(u) - ∂x(w))   # vorticity in the xz plane

simulation.output_writers[:fields] = NetCDFWriter(model,
    (; u, w, ω),
    schedule           = TimeInterval(0.5),
    filename           = "hill_flow.nc",
    overwrite_existing = true)

run!(simulation)
