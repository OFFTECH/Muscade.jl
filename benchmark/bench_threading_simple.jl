using Pkg
Pkg.activate(".")
Pkg.instantiate()

using Muscade
using BenchmarkTools
using LinearAlgebra

println("=== Threading Benchmark for SweepX Assembly ===")
println()

function run_benchmark(nthreads::Int)
    println("Running with JULIA_NUM_THREADS=$nthreads")
    
    model = Model(:TestModel)
    node1 = addnode!(model, 𝕣[0,0,0])
    node2 = addnode!(model, 𝕣[1,0,0])
    
    @functor with() sea(t,x) = SVector(1.,0.)*t
    @functor with() sky(t,x) = SVector(0.,10.)
    
    addelement!(model, QuickFix, [node1]; field=:x, value=0.)
    addelement!(model, QuickFix, [node1]; field=:y, value=0.)
    addelement!(model, QuickFix, [node2]; field=:x, value=sea)
    addelement!(model, QuickFix, [node2]; field=:y, value=sky)
    
    initialstate = initialize!(model)
    
    b = @benchmark solve(SweepX{2}; initialstate=deepcopy(initialstate), time=0:0.1:1., verbose=false) samples=10 evals=1
    
    return b
end

println("Note: This benchmark uses the test model from TestSweepX0.jl")
println("To compare threading, run with different JULIA_NUM_THREADS values")
println()

b = run_benchmark(Threads.nthreads())
display(b)
println()
