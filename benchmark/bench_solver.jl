using Pkg
Pkg.activate(dirname(@__DIR__))
Pkg.instantiate()

using Muscade
using Muscade.Toolbox
using Muscade: NodID, Hold, DofLoad
using BenchmarkTools
using Printf
using StaticArrays

println("=" ^ 60)
println("Benchmark: Solver memory allocation")
println("=" ^ 60)
println()

# Define constant load
@functor with() load(t) = -1000.0

# Create a simple beam model for benchmarking
function create_beam_model(nelements::Int)
    model = Model(:BenchmarkBeam)
    
    L = 10.0  # Total length
    dx = L / nelements
    
    # Create nodes - addnode! returns NodID
    nodes = Vector{NodID}()
    for i = 1:nelements+1
        x = (i-1) * dx
        push!(nodes, addnode!(model, [x, 0.0, 0.0]))
    end
    
    # Create beam elements
    mat = BeamCrossSection(EA=1e8, EI₂=1e6, EI₃=1e6, GJ=1e5, μ=100.0, ι₁=10.0)
    for i = 1:nelements
        addelement!(model, EulerBeam3D, [nodes[i], nodes[i+1]]; mat=mat)
    end
    
    # Boundary conditions - clamp at first node
    for field ∈ [:t1, :t2, :t3, :r1, :r2, :r3]
        addelement!(model, Hold, [nodes[1]]; field=field)
    end
    
    # Load at tip
    addelement!(model, DofLoad, [nodes[end]]; field=:t3, value=load)
    
    return model
end

println("Creating model with 100 beam elements...")
model = create_beam_model(100)
initialstate = initialize!(model)

println("Running SweepX{0} solver benchmark...")
println()

# Warm up
states = solve(SweepX{0}; initialstate=initialstate, time=[0.0], verbose=false)

# Benchmark
b = @benchmark solve(SweepX{0}; initialstate=$initialstate, time=[0.0], verbose=false) samples=10 evals=1
display(b)
println()

# Memory allocation test
println("Memory allocation test:")
alloc = @allocated solve(SweepX{0}; initialstate=initialstate, time=[0.0], verbose=false)
println("  Total allocations: $alloc bytes")
println()

println("=" ^ 60)
println("Benchmark complete")
println("=" ^ 60)
