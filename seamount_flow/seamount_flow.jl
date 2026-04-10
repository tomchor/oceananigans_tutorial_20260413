using Oceananigans
using Oceananigans.Solvers: ConjugateGradientPoissonSolver, fft_poisson_solver
using NCDatasets
using Printf

# =============================================================================
# Flow past a Gaussian seamount
#
# A 2-D (xz) nonhydrostatic simulation with open east/west boundaries.
# A background flow U∞ enters from the west, goes over a seamount, and exits
# to the east.  The immersed boundary method handles the topography.
#
# Stratification (N² > 0) can be added to show internal wave generation.
# =============================================================================

# --- Physical parameters ---
Lx = 20.0       # m, domain length
H  = 2.0        # m, domain depth
U∞ = 1.0        # m/s, inflow velocity

# Seamount geometry
x₀ = Lx / 3    # center position
h₀ = 0.6H      # peak height above the bottom  (= 0.6 * H above z = -H)
σ  = Lx / 10   # horizontal half-width

seamount(x, y, p) = p.h₀ * exp(-((x - p.x₀) / p.σ)^2) - p.H   # returns z_bottom(x)

# --- Grid ---
Nx, Nz = 128, 32

underlying_grid = RectilinearGrid(size     = (Nx, Ny, Nz),
                                  x        = (-Lx/2, Lx/2),
                                  y        = (-Ly/2, Lx/2),
                                  z        = (0, H),
                                  topology = (Bounded, Flat, Bounded),
                                  halo     = (6, 6))

seamount_params = (; x₀, h₀, σ, H)
grid = ImmersedBoundaryGrid(underlying_grid, GridFittedBottom((x, y) -> seamount(x, y, seamount_params)))

# --- Pressure solver (required for non-periodic x topology) ---
pressure_solver = ConjugateGradientPoissonSolver(grid;
    maxiter      = 10,
    preconditioner = fft_poisson_solver(underlying_grid))

# --- Open boundary conditions on east/west faces ---
u_bcs = FieldBoundaryConditions(west = OpenBoundaryCondition(U∞),
                                east = OpenBoundaryCondition(U∞))

# --- Model ---
model = NonhydrostaticModel(grid;
                            pressure_solver,
                            boundary_conditions = (u=u_bcs,),
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
    @info @sprintf("t = %s, Δt = %s, max|u| = %.3f, wall time = %s",
                   prettytime(time(sim)), prettytime(sim.Δt), maximum(abs, u), elapsed)
    wall_clock[] = time_ns()
end
add_callback!(simulation, progress, IterationInterval(100))

# --- Output ---
u, v, w = model.velocities
ω = Field(∂z(u) - ∂x(w))    # vorticity in the xz plane

simulation.output_writers[:fields] = NetCDFWriter(model,
    (; u, w, ω),
    schedule           = TimeInterval(0.5),
    filename           = "seamount_flow.nc",
    overwrite_existing = true)

run!(simulation)
