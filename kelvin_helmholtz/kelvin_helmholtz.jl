using Oceananigans
using NCDatasets
using Printf
using Oceanostics.FlowDiagnostics: StrainRateTensorModulus

# =============================================================================
# Kelvin-Helmholtz instability (2D, xz plane) — implicit LES
#
# A shear flow u = U·tanh(z/h) with stable stratification b = B₀·tanh(z/h).
# WENO(order=5) dissipation replaces explicit viscosity (implicit LES).
#
# Increase Nz for better-resolved billows.
# =============================================================================

# --- Physical parameters ---
U                      = 1.0    # velocity profile amplitude (m/s)
Ri                     = 0.1    # Richardson number (must be < 0.25 for instability)
h                      = 1.0    # shear/buoyancy layer half-width (m)
perturbation_amplitude = 0.05   # perturbation amplitude

B₀ = U^2 * Ri / h   # buoyancy amplitude

# --- Most unstable KH wavenumber (Michalke 1964) ---
k_max = 0.4446 / h
λ_max = 2π / k_max

# --- Domain ---
Lx = λ_max
Lz = 15 * h

# --- Grid ---
Nz = 256
Nx = round(Int, Nz * (Lx / Lz) / 2)   # cell aspect ratio Δx/Δz ≈ 2

@info @sprintf("Most unstable KH wavenumber: k_max = %.4f  (λ_max = %.2f, Lx = %.1f)", k_max, λ_max, Lx)

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
uᵢ(x, z) = U * tanh(z / h)
bᵢ(x, z) = B₀ * tanh(z / h)
wᵢ(x, z) = perturbation_amplitude * abs(randn()) * exp(-z^2) * sin(x * k_max - π) # Nice and centered eye
set!(model, u=uᵢ, b=bᵢ, w=wᵢ)

# --- Simulation ---
Δt = 0.2 * minimum_xspacing(grid) / maximum(model.velocities.u)
simulation = Simulation(model; Δt=Δt, stop_time=200)
conjure_time_step_wizard!(simulation, cfl=0.8, IterationInterval(5))

function progress(sim)
    w = sim.model.velocities.w
    percent = 100 * time(sim) / sim.stop_time
    @info @sprintf("t = %.4f, Δt = %.4f, max|w| = %.3f (%.1f%% complete)", time(sim), sim.Δt, maximum(abs, w), percent)
end
add_callback!(simulation, progress, IterationInterval(100))

# --- Output ---
u, v, w = model.velocities
b = model.tracers.b
ω = Field(∂z(u) - ∂x(w))   # vorticity in the xz plane
S = StrainRateTensorModulus(model)

simulation.output_writers[:fields] = NetCDFWriter(model,
    (; u, w, b, ω, S),
    schedule           = TimeInterval(2.0),
    filename           = "kelvin_helmholtz.nc",
    overwrite_existing = true)

run!(simulation)
