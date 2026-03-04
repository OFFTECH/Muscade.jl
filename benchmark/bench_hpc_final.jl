using Pkg
Pkg.activate(".")

using Muscade, BenchmarkTools

include("../test/SomeElements.jl")

function run_benchmark(NELEMENTS)
    model = Model(:BenchModel)
    nodes = [addnode!(model, 𝕣[]) for _ in 1:NELEMENTS]
    for i in 1:NELEMENTS
        addelement!(model, SdofOscillator, [nodes[i]]; K₁=1., K₂=0.3, C₁=0.1, C₂=0., M₁=1., M₂=0.)
    end

    mystate = initialize!(model; time=0.)
    for i in 1:NELEMENTS
        mystate = setdof!(mystate, [1.0]; field=:tx1, nodID=[nodes[i]], order=0)
    end

    solve(SweepX{2}; initialstate=deepcopy(mystate), time=0.01:0.1:0.11, verbose=false)

    b = @benchmark solve(SweepX{2}; initialstate=deepcopy($mystate), time=0.01:0.1:1.0, verbose=false) samples=10 evals=1
    
    return b
end

println("=" ^ 60)
println("HPC Benchmark Results")
println("=" ^ 60)
println("Threads: $(Threads.nthreads())")
println()

for n in [200, 500, 1000]
    println("-" ^ 40)
    println("Elements: $n")
    b = run_benchmark(n)
    println("  Min:    $(round(minimum(b).time/1e6, digits=2)) ms")
    println("  Median: $(round(median(b).time/1e6, digits=2)) ms")
    println("  Mean:   $(round(mean(b).time/1e6, digits=2)) ms")
    println("  Memory: $(round(BenchmarkTools.memory(b)/1024/1024, digits=2)) MiB")
    println()
end
