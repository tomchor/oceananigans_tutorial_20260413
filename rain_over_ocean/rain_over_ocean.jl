using Breeze
using Oceananigans.Units

using AtmosphericProfilesLibrary
using CloudMicrophysics
using Printf
using Random

# =============================================================================
# Precipitating shallow cumulus convection (RICO, 2D xz) — Breeze / Oceananigans
#
# Rain in Cumulus over the Ocean intercomparison case (vanZanten et al. 2011).
# Trade-wind shallow cumulus with warm-rain microphysics, large-scale subsidence,
# prescribed radiative cooling and moisture tendencies, and bulk surface fluxes.
# =============================================================================

Random.seed!(42)

# --- Grid ---
Nx = 128
Nz = 100

grid = RectilinearGrid(; size = (Nx, Nz),
                         x        = (0, 12800),   # 12.8 km
                         z        = (0, 4000),    # 4 km
                         halo     = (5, 5),
                         topology = (Periodic, Flat, Bounded))

# --- Reference state ---
constants       = ThermodynamicConstants()
reference_state = ReferenceState(grid, constants;
                                 surface_pressure      = 101540,
                                 potential_temperature = 297.9)
dynamics = AnelasticDynamics(reference_state)

# --- Bulk surface fluxes (vanZanten et al. 2011, eqs. 1–4) ---
Cᵀ = 1.094e-3   # sensible heat transfer coefficient
Cᵛ = 1.133e-3   # moisture transfer coefficient
T₀ = 299.8      # sea-surface temperature (K)

ρe_bcs  = FieldBoundaryConditions(bottom = BulkSensibleHeatFlux(coefficient = Cᵀ, surface_temperature = T₀))
ρqᵉ_bcs = FieldBoundaryConditions(bottom = BulkVaporFlux(coefficient = Cᵛ, surface_temperature = T₀))

# --- Sponge layer (damps w in upper 500 m) ---
# GaussianMask only has a 3-arg method (x,y,z) so we use a plain 2-arg function for 2D
sponge_mask(x, z) = exp(-((z - 3500)^2) / (2 * 500^2))
sponge = Relaxation(rate = 1/8, mask = sponge_mask)

# --- Large-scale subsidence ---
FT = eltype(grid)
wˢ = Field{Nothing, Nothing, Face}(grid)
set!(wˢ, z -> AtmosphericProfilesLibrary.Rico_subsidence(FT)(z))
subsidence = SubsidenceForcing(wˢ)

# --- Geostrophic wind forcing (f = 4.5×10⁻⁵ s⁻¹, ≈18°N) ---
coriolis  = FPlane(f = 4.5e-5)
geostrophic = geostrophic_forcings(z -> AtmosphericProfilesLibrary.Rico_geostrophic_ug(FT)(z),
                                   z -> AtmosphericProfilesLibrary.Rico_geostrophic_vg(FT)(z))

# --- Large-scale moisture tendency ---
ρᵣ = reference_state.density
∂t_ρqᵉ = Field{Nothing, Nothing, Center}(grid)
set!(∂t_ρqᵉ, z -> AtmosphericProfilesLibrary.Rico_dqtdt(FT)(z))
set!(∂t_ρqᵉ, ρᵣ * ∂t_ρqᵉ)
∂t_ρqᵉ_forcing = Forcing(∂t_ρqᵉ)

# --- Radiative cooling (−2.5 K/day, uniform) ---
∂t_ρθ = Field{Nothing, Nothing, Center}(grid)
set!(∂t_ρθ, ρᵣ * (-2.5 / day))
ρθ_forcing = Forcing(∂t_ρθ)

# --- Assemble forcing and boundary conditions ---
forcing = (ρu  = (subsidence, geostrophic.ρu),
           ρv  = (subsidence, geostrophic.ρv),
           ρw  = sponge,
           ρqᵉ = (subsidence, ∂t_ρqᵉ_forcing),
           ρθ  = (subsidence, ρθ_forcing))

boundary_conditions = (ρe = ρe_bcs, ρqᵉ = ρqᵉ_bcs)

# --- One-moment cloud microphysics (autoconversion + accretion) ---
BreezeCloudMicrophysicsExt = Base.get_extension(Breeze, :BreezeCloudMicrophysicsExt)
using .BreezeCloudMicrophysicsExt: OneMomentCloudMicrophysics

microphysics    = OneMomentCloudMicrophysics(cloud_formation = SaturationAdjustment(equilibrium = WarmPhaseEquilibrium()))
weno            = WENO(order = 5)
bpweno          = WENO(order = 5, bounds = (0, 1))   # bounds-preserving for moisture
scalar_advection = (ρθ = weno, ρqᵉ = bpweno, ρqᶜˡ = bpweno, ρqʳ = bpweno)

model = AtmosphereModel(grid; dynamics, coriolis, microphysics,
                        momentum_advection = weno,
                        scalar_advection,
                        forcing, boundary_conditions)

# --- Initial conditions ---
θˡⁱ₀ = AtmosphericProfilesLibrary.Rico_θ_liq_ice(FT)
qᵗ₀  = AtmosphericProfilesLibrary.Rico_q_tot(FT)
u₀   = AtmosphericProfilesLibrary.Rico_u(FT)
v₀   = AtmosphericProfilesLibrary.Rico_v(FT)

zϵ = 1500   # add small random perturbations only below 1500 m to seed convection
θᵢ(x, z) = θˡⁱ₀(z) + 1e-2 * (rand() - 0.5) * (z < zϵ)
qᵢ(x, z) = qᵗ₀(z)
uᵢ(x, z) = u₀(z)
vᵢ(x, z) = v₀(z)

set!(model, θ = θᵢ, qᵗ = qᵢ, u = uᵢ, v = vᵢ)

# --- Simulation ---
simulation = Simulation(model; Δt = 2, stop_time = 5hour)
conjure_time_step_wizard!(simulation, cfl = 0.7)

θ   = liquid_ice_potential_temperature(model)
qᶜˡ = model.microphysical_fields.qᶜˡ
qʳ  = model.microphysical_fields.qʳ

wall_clock = Ref(time_ns())

function progress(sim)
    wmax    = maximum(abs, sim.model.velocities.w)
    elapsed = 1e-9 * (time_ns() - wall_clock[])
    percent = 100 * time(sim) / sim.stop_time
    @info @sprintf("t = %s, Δt = %s, wall: %s, max|w| = %.2e m/s, max(qᶜˡ) = %.2e, max(qʳ) = %.2e (%.1f%% done)",
                   prettytime(sim), prettytime(sim.Δt), prettytime(elapsed),
                   wmax, maximum(qᶜˡ), maximum(qʳ), percent)
end

add_callback!(simulation, progress, IterationInterval(100))

# --- Output ---
u, v, w = model.velocities

simulation.output_writers[:fields] = JLD2Writer(model,
    (; u, w, θ, qᶜˡ, qʳ),
    schedule           = TimeInterval(2minutes),
    filename           = "rain_over_ocean.jld2",
    overwrite_existing = true)

run!(simulation)
