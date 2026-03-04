using Pkg
Pkg.activate(".")

using Muscade, BenchmarkTools

include("../test/SomeElements.jl")

function run_benchmark()
    NELEMENTS = 1000

    model = Model(:BenchModel)
    nodes = [addnode!(model, 𝕣[]) for _ in 1:NELEMENTS]
    for i in 1:NELEMENTS
        addelement!(model, SdofOscillator, [nodes[i]]; K₁=1., K₂=0.3, C₁=0.1, C₂=0., M₁=1., M₂=0.)
    end

    mystate = initialize!(model; time=0.)
    for i in 1:NELEMENTS
        mystate = setdof!(mystate, [1.0]; field=:tx1, nodID=[nodes[i]], order=0)
    end

    println("Model: $NELEMENTS SdofOscillator elements")
    println("Threading threshold: 50")
    println("Threads: $(Threads.nthreads())")
    println()

    # Warmup - time starts at 0.01, not 0
    solve(SweepX{2}; initialstate=deepcopy(mystate), time=0.01:0.1:0.11, verbose=false)

    # Benchmark
    b = @benchmark solve(SweepX{2}; initialstate=deepcopy($mystate), time=0.01:0.1:1.0, verbose=false) samples=5 evals=1
    display(b)
    println()
end

run_benchmark()
