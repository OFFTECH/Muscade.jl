using Pkg
Pkg.activate(".")
Pkg.instantiate()

using Muscade
using Muscade.Toolbox: Bar3D, AxisymmetricBarCrossSection
using BenchmarkTools
using LinearAlgebra
using Random

Random.seed!(42)

@functor with() load100(t) = 1e6 * sin(2π * t)
@functor with() load1000(t) = 1e6 * sin(2π * t)

function create_truss_model_100()
    nelements = 100
    model = Model()
    
    nodes = [addnode!(model, [Float64(i-1), 0., 0.]) for i in 1:nelements+1]
    
    mat = AxisymmetricBarCrossSection(EA=2e9, μ=7850.)
    
    for i in 1:nelements
        addelement!(model, Bar3D, [nodes[i], nodes[i+1]]; mat=mat)
    end
    
    addelement!(model, Hold, [nodes[1]], field = :x)
    addelement!(model, Hold, [nodes[1]], field = :y)
    addelement!(model, Hold, [nodes[1]], field = :z)
    
    addelement!(model, DofLoad, [nodes[end]]; field = :y, value = load100)
    
    return model
end

function create_truss_model_1000()
    nelements = 1000
    model = Model()
    
    nodes = [addnode!(model, [Float64(i-1), 0., 0.]) for i in 1:nelements+1]
    
    mat = AxisymmetricBarCrossSection(EA=2e9, μ=7850.)
    
    for i in 1:nelements
        addelement!(model, Bar3D, [nodes[i], nodes[i+1]]; mat=mat)
    end
    
    addelement!(model, Hold, [nodes[1]], field = :x)
    addelement!(model, Hold, [nodes[1]], field = :y)
    addelement!(model, Hold, [nodes[1]], field = :z)
    
    addelement!(model, DofLoad, [nodes[end]]; field = :y, value = load1000)
    
    return model
end

println("Creating models...")
model_100 = create_truss_model_100()
model_1000 = create_truss_model_1000()

println("\nInitializing...")
state_100 = initialize!(model_100)
state_1000 = initialize!(model_1000)

println("\n=== Benchmark: 100 elements (threading threshold = 50) ===")
println("Threads: $(Threads.nthreads())")

b100 = @benchmark solve(SweepX{2}; initialstate=state_100, time=0:0.01:0.1, verbose=false) samples=5 evals=1
println("\n100 elements:")
display(b100)

println("\n\n=== Benchmark: 1000 elements (threading threshold = 50) ===")
b1000 = @benchmark solve(SweepX{2}; initialstate=state_1000, time=0:0.01:0.1, verbose=false) samples=3 evals=1
println("\n1000 elements:")
display(b1000)

println("\n\n=== Speedup analysis ===")
println("Threads available: $(Threads.nthreads())")
println("For 100 elements: threading $(100 >= 50 ? "ENABLED" : "disabled")")
println("For 1000 elements: threading $(1000 >= 50 ? "ENABLED" : "disabled")")
