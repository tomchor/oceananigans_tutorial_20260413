using Oceananigans
using Oceananigans.Units
using NCDatasets
using Printf

# =============================================================================
# Atmospheric free convection
#
# A weakly stratified atmosphere heated from below. No mean pressure gradient
# — flow is driven purely by buoyancy (free convection).
# Quadratic bulk drag at the surface.
#
# Temperature (potential temperature θ) is the active tracer, coupled to
# buoyancy via SeawaterBuoyancy with an ideal-gas linear EOS (β = 0).
# =============================================================================

# --- Physical parameters ---
θ₀  = 300.0    # K, reference potential temperature
g   = 9.81     # m/s², gravitational acceleration
α   = 1 / θ₀   # K⁻¹, thermal expansion coefficient (ideal gas: α = 1/θ₀)
N²  = 1e-4     # s⁻², initial buoyancy frequency (weak stratification)
Qᵀ  = 5e-2    # K m/s, kinematic surface heat flux (positive = warming)
z₀  = 0.1meters  # aerodynamic roughness length (moderate terrain / farmland)

# --- Domain ---
Lx = Ly = 4kilometers
H  = 2kilometers

# --- Grid ---
Nx = Ny = 64
Nz = 32

grid = RectilinearGrid(size     = (Nx, Ny, Nz),
                       x        = (0, Lx),
                       y        = (0, Ly),
                       z        = (0, H),
                       topology = (Periodic, Periodic, Bounded))

# --- Drag coefficient from law of the wall ---
const κᵛᵏ = 0.4   # von Kármán constant
z₁ = minimum_zspacing(grid, Center(), Center(), Center()) / 2
Cd = (κᵛᵏ / log(z₁ / z₀))^2
@info "z₁ = $(z₁) m,  Cd = $(round(Cd, sigdigits=3))"

# --- Boundary conditions ---
drag  = BulkDrag(coefficient=Cd)
u_bcs = FieldBoundaryConditions(bottom = drag)
v_bcs = FieldBoundaryConditions(bottom = drag)

# Surface heat flux at z = 0: positive flux warms the surface layer
T_bcs = FieldBoundaryConditions(bottom = FluxBoundaryCondition(Qᵀ))

# --- Buoyancy: temperature via SeawaterBuoyancy with ideal-gas linear EOS ---
# b = g · α · (T - T_ref),  β = 0 (dry atmosphere, no salinity)
buoyancy = SeawaterBuoyancy(equation_of_state = LinearEquationOfState(thermal_expansion = α,
                                                                       haline_contraction = 0),
                             constant_salinity = 0)

# --- Model ---
model = NonhydrostaticModel(grid;
                            closure             = SmagorinskyLilly(),
                            advection           = WENO(order=5),
                            timestepper         = :RungeKutta3,
                            buoyancy            = buoyancy,
                            tracers             = :T,
                            boundary_conditions = (u=u_bcs, v=v_bcs, T=T_bcs))

# --- Initial conditions ---
# Weak linear stratification + small-amplitude noise to seed convection
dθdz = N² / (g * α)   # K/m, consistent with chosen N²
θᵢ(x, y, z) = θ₀ + dθdz * z + 0.01 * randn()
set!(model, T=θᵢ)

# --- Simulation ---
w★  = (g / θ₀ * Qᵀ * H)^(1/3)   # Deardorff velocity scale (m/s)
Δt₀ = 0.1 * minimum_zspacing(grid) / w★
simulation = Simulation(model; Δt=Δt₀, stop_time=4hours)
conjure_time_step_wizard!(simulation, cfl=0.8, IterationInterval(5))

wall_clock = Ref(time_ns())
function progress(sim)
    w = sim.model.velocities.w
    T = sim.model.tracers.T
    elapsed = prettytime(1e-9 * (time_ns() - wall_clock[]))
    @info @sprintf("t = %s, Δt = %s, max|w| = %.2f m/s, T_sfc = %.3f K, wall time = %s",
                   prettytime(time(sim)), prettytime(sim.Δt), maximum(abs, w),
                   maximum(interior(T, :, :, 1)), elapsed)
    wall_clock[] = time_ns()
end
add_callback!(simulation, progress, IterationInterval(100))

# --- Output ---
u, v, w = model.velocities
T = model.tracers.T

simulation.output_writers[:midlevel] = NetCDFWriter(model,
    (; u, v, w, T),
    schedule           = TimeInterval(2minutes),
    filename           = "free_convection.nc",
    indices            = (:, :, Nz÷2),     # mid-level (z ≈ H/2) horizontal slice
    overwrite_existing = true)

simulation.output_writers[:xz_slice] = NetCDFWriter(model,
    (; u, v, w, T),
    schedule           = TimeInterval(2minutes),
    filename           = "free_convection_xz.nc",
    indices            = (:, Ny÷2, :),     # xz slice at mid-domain y
    overwrite_existing = true)

run!(simulation)
