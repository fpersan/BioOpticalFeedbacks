using Oceananigans
using Oceananigans.Units
using Oceananigans.TurbulenceClosures: CATKEVerticalDiffusivity
using GLMakie

grid = RectilinearGrid(size=64, z=(-256, 0), topology=(Flat, Flat, Bounded))
coriolis = FPlane(f=1e-4)
closure = CATKEVerticalDiffusivity()

# Initial buoyancy frequency
N² = 1e-6

# Surface buoyancy flux (convection is positive flux upward)
Qᵇ = +1e-8

# Surface momentum flux (convection is positive flux upward)
Qᵘ = -2e-4 #

b_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(Qᵇ))
u_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(Qᵘ))

model = HydrostaticFreeSurfaceModel(; grid, closure, coriolis,
                                    tracers = (:b, :e),
                                    buoyancy = BuoyancyTracer(),
                                    boundary_conditions = (; b=b_bcs, u=u_bcs))

bᵢ(z) = N² * z
set!(model, b=bᵢ, e=1e-6)

simulation = Simulation(model, Δt=10minutes, stop_iteration=1000)

closurename = string(nameof(typeof(closure)))

diffusivities = (κᵘ = model.diffusivity_fields.κᵘ,
                 κᶜ = model.diffusivity_fields.κᶜ)

outputs = merge(model.velocities, model.tracers, diffusivities)

simulation.output_writers[:fields] = JLD2OutputWriter(model, outputs,
                                                      schedule = TimeInterval(10minutes),
                                                      filename = "windy_convection.jld2",
                                                      overwrite_existing = true)

progress(sim) = @info string("Iter: ", iteration(sim), " t: ", prettytime(sim))
simulation.callbacks[:progress] = Callback(progress, IterationInterval(100))

@info "Running a simulation of $model..."

run!(simulation)

ut = FieldTimeSeries("windy_convection.jld2", "u")
bt = FieldTimeSeries("windy_convection.jld2", "b")
times = bt.times
Nt = length(times)

fig = Figure()
axu = Axis(fig[1, 1])
axb = Axis(fig[1, 2])

slider = Slider(fig[2, 1:2], range=1:Nt, startvalue=1)
n = slider.value

bn = @lift interior(bt[$n], 1, 1, :)
un = @lift interior(ut[$n], 1, 1, :)

z = znodes(b)

lines!(axu, un, z)
lines!(axb, bn, z)

display(fig)

