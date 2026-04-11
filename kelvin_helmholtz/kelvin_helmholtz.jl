using Oceananigans
using NCDatasets
using Printf

# =============================================================================
# Kelvin-Helmholtz instability (2D, xz plane)
#
# A shear flow u = tanh(z) with stable stratification b = h·Ri·tanh(z/h).
# A small perturbation at the most unstable wavenumber triggers the instability.
#
# Increase Nz for better-resolved billows.
# =============================================================================

# --- Physical parameters ---
Ri  = 0.1     # Richardson number (must be < 0.25 for instability)
h   = 0.25    # buoyancy layer thickness
perturbation_amplitude = 0.01

# --- Grid ---
Lx, Lz = 10.0, 14.0
Nz = 256
Nx = round(Int, Nz * Lx / Lz / 2)   # cell aspect ratio Δx/Δz ≈ 2

grid = RectilinearGrid(size     = (Nx, Nz),
                       x        = (-Lx/2, Lx/2),
                       z        = (-Lz/2, Lz/2),
                       topology = (Periodic, Flat, Bounded))

# --- Model (implicit LES: WENO dissipation replaces explicit viscosity) ---
model = NonhydrostaticModel(grid;
                            advection = WENO(order=5),
                            buoyancy  = BuoyancyTracer(),
                            tracers   = :b)

# --- Initial conditions ---
# Most unstable KH wavenumber (Michalke 1964; Hazel 1972 stratification correction)
k_max = 0.4446 * sqrt(max(0.0, 1 - 4*Ri))
λ_max = 2π / k_max
@info @sprintf("Most unstable KH wavenumber: k_max = %.4f  (λ_max = %.2f, Lx = %.1f)",
               k_max, λ_max, Lx)

uᵢ(x, z) = tanh(z) + perturbation_amplitude * sin(2π/Lx * x) * exp(-z^2 / 2)
bᵢ(x, z) = h * Ri * tanh(z / h)
wᵢ(x, z) = perturbation_amplitude * cos(2π/Lx * x) * exp(-z^2 / 2)
set!(model, u=uᵢ, b=bᵢ, w=wᵢ)

# --- Simulation ---
simulation = Simulation(model; Δt=0.01, stop_time=200.0)
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
b = model.tracers.b
ω = Field(∂z(u) - ∂x(w))   # vorticity in the xz plane

simulation.output_writers[:fields] = NetCDFWriter(model,
    (; u, w, b, ω),
    schedule           = TimeInterval(2.0),
    filename           = "kelvin_helmholtz.nc",
    overwrite_existing = true)

run!(simulation)
