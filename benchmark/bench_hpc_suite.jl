# HPC Benchmark Suite for# This file provides standardized benchmarks for# Muscade.jl HPC performance testing

using Pkg
Pkg.activate(".")

using Muscade
using BenchmarkTools
using Printf
using Random

# ============================================================
# Benchmark Runner
# ============================================================

const RESULTS = Dict{String, Any}()

function print_header(title)
    println("\n" * "=" ^ 50)
    println("  $title")
    println("  " * "=" ^ 50)
end

function print_result(name, time_ms, memory_mb)
    @printf("  %-25s: %8.2f ms | %6.2f MB\n", name, time_ms, memory_mb)
end

# ============================================================
# Benchmark 1: Assembly Operations
# ============================================================

function bench_assembly_zero!()
    n = 1000
    v = randn(n)
    sparse_v = sprand(n, 0.01)
    
    b = @benchmark zero!(v) samples=100
    print_result("zero! (dense)", minimum(b).time / 1e6, BenchmarkTools.memory(b) / 1024 / 1024)
    
    b = @benchmark zero!(sparse_v) samples=100
    print_result("zero! (sparse)", minimum(b).time / 1e6, BenchmarkTools.memory(b) / 1024 / 1024)
end

# ============================================================
# Benchmark 2: DofGroup Operations
# ============================================================

function bench_dofgroup()
    n = 500
    model = Model(:BenchModel)
    
    # Create nodes
    nodes = [addnode!(model, 𝕣[Float64(i)]) for i in 1:n]
    
    # Create simple elements
    include("../test/SomeElements.jl")
    for i in 1:n
        addelement!(model, SdofOscillator, [nodes[i]]; K₁=1., K₂=0.3, C₁=0.1, C₂=0., M₁=1., M₂=0.)
    end
    
    state = initialize!(model; time=0.)
    
    # Benchmark getdof!
    b = @benchmark getdof(state; field=:tx1, nodID=[nodes[1]]) samples=100
    print_result("getdof", minimum(b).time / 1e6, BenchmarkTools.memory(b) / 1024 / 1024)
end

# ============================================================
# Benchmark 3: Solver Performance
# ============================================================

function bench_solver(nelements::Int)
    model = Model(:BenchModel)
    nodes = [addnode!(model, 𝕣[]) for _ in 1:nelements]
    
    include("../test/SomeElements.jl")
    for i in 1:nelements
        addelement!(model, SdofOscillator, [nodes[i]]; K₁=1., K₂=0.3, C₁=0.1, C₂=0., M₁=1., M₂=0.)
    end
    
    state = initialize!(model; time=0.)
    for i in 1:nelements
        state = setdof!(state, [1.0]; field=:tx1, nodID=[nodes[i]], order=0)
    end
    
    # Warmup
    solve(SweepX{2}; initialstate=deepcopy(state), time=0.01:0.1, verbose=false)
    
    # Benchmark
    b = @benchmark solve(SweepX{2}; initialstate=deepcopy($state), time=0.01:0.1, verbose=false) samples=10 evals=1
    
    return (nelements, minimum(b).time / 1e6, median(b).time / 1e6, BenchmarkTools.memory(b) / 1024 / 1024)
end

# ============================================================
# Main
# ============================================================

println("=" ^ 60)
println("MUSCADE.jl HPC BENCHMARK SUITE")
println("=" ^ 60)
println("Threads: $(Threads.nthreads())")
println("Julia version: $(VERSION)")

print_header("1. Assembly Operations")
bench_assembly_zero!()

print_header("2. DofGroup Operations")
bench_dofgroup()

print_header("3. Solver Scaling")
println("  Elements | Min (ms) | Median (ms) | Memory (MB)")
println("  " * "-" ^ 40)

for n in [100, 500, 1000]
    nelems, t_min, t_med, mem = bench_solver(n)
    @printf("  %-7d | %7.1f | %7.1f | %6.2f\n", n, t_min / 1e6, t_med / 1e6, mem)
end

println("\n" * "=" ^ 60)
println("BENCHMARK COMPLETE")
println("  " * "=" ^ 60)
