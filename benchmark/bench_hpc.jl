# HPC Performance Benchmark Suite for Muscade.jl
# This benchmark measures key performance indicators identified in the HPC analysis

using Pkg
Pkg.activate(dirname(@__DIR__))
Pkg.instantiate()

using Muscade
using SparseArrays
using LinearAlgebra
using BenchmarkTools
using Printf

println("=" ^ 70)
println("Muscade.jl HPC Performance Benchmark Suite")
println("=" ^ 70)
println()

# Helper function to format benchmark results
function format_benchmark(name, b)
    mean_time = mean(b).time / 1e6  # Convert to ms
    std_time = std(b).time / 1e6
    min_time = minimum(b).time / 1e6
    memory = b.memory
    allocs = b.allocs
    @printf "  %-35s: %8.3f ms ± %6.3f ms (min: %8.3f ms, allocs: %d, memory: %s)\n" name mean_time std_time min_time allocs format_bytes(memory)
    return (name=name, mean=mean_time, std=std_time, min=min_time, allocs=allocs, memory=memory)
end

function format_bytes(bytes)
    if bytes < 1024
        return @sprintf("%d B", bytes)
    elseif bytes < 1024^2
        return @sprintf("%.1f KB", bytes / 1024)
    elseif bytes < 1024^3
        return @sprintf("%.1f MB", bytes / 1024^2)
    else
        return @sprintf("%.1f GB", bytes / 1024^3)
    end
end

results = []

# ============================================================================
# Section 1: zero! function benchmarks
# ============================================================================
println("=" ^ 70)
println("Section 1: zero! function performance")
println("=" ^ 70)
println()

println("1.1 Dense Vector zero! (n=10_000)")
v1 = zeros(10_000)
b1 = @benchmark Muscade.zero!($v1) samples=1000 evals=10
push!(results, format_benchmark("zero! dense 10k", b1))

println("1.2 Dense Vector zero! (n=100_000)")
v2 = zeros(100_000)
b2 = @benchmark Muscade.zero!($v2) samples=1000 evals=10
push!(results, format_benchmark("zero! dense 100k", b2))

println("1.3 Sparse Matrix zero! (1000x1000, ~5000 nnz)")
sp1 = sprand(1000, 1000, 0.005)
b3 = @benchmark Muscade.zero!($sp1) samples=1000 evals=10
push!(results, format_benchmark("zero! sparse 5k nnz", b3))

println("1.4 Sparse Matrix zero! (10000x10000, ~50000 nnz)")
sp2 = sprand(10000, 10000, 0.0005)
b4 = @benchmark Muscade.zero!($sp2) samples=1000 evals=10
push!(results, format_benchmark("zero! sparse 50k nnz", b4))

println()

# ============================================================================
# Section 2: Linear algebra benchmarks
# ============================================================================
println("=" ^ 70)
println("Section 2: Linear algebra performance")
println("=" ^ 70)
println()

println("2.1 Sparse LU factorization (1000x1000)")
A1 = sprand(1000, 1000, 0.01) + sparse(I, 1000, 1000) * 100  # Ensure non-singular
b_lu1 = @benchmark lu($A1) samples=100 evals=1
push!(results, format_benchmark("LU factorization 1k", b_lu1))

println("2.2 Sparse LU factorization (5000x5000)")
A2 = sprand(5000, 5000, 0.002) + sparse(I, 5000, 5000) * 100
b_lu2 = @benchmark lu($A2) samples=50 evals=1
push!(results, format_benchmark("LU factorization 5k", b_lu2))

println("2.3 Sparse solve with existing factorization (1000x1000)")
F1 = lu(A1)
rhs1 = rand(1000)
b_solve1 = @benchmark ldiv!($rhs1, $F1, copy($rhs1)) samples=1000 evals=10
push!(results, format_benchmark("ldiv! with LU 1k", b_solve1))

println("2.4 Sparse solve with allocation (1000x1000) - baseline")
b_solve_alloc = @benchmark $F1 \ $rhs1 samples=1000 evals=10
push!(results, format_benchmark("F\\b allocation 1k", b_solve_alloc))

println()

# ============================================================================
# Section 3: Memory allocation benchmarks
# ============================================================================
println("=" ^ 70)
println("Section 3: Memory allocation tracking")
println("=" ^ 70)
println()

println("3.1 zero! allocations")
v_alloc = zeros(100_000)
alloc1 = @allocated Muscade.zero!(v_alloc)
println("  Allocations for zero! on 100k vector: $(format_bytes(alloc1))")

sp_alloc = sprand(10000, 10000, 0.0005)
alloc2 = @allocated Muscade.zero!(sp_alloc)
println("  Allocations for zero! on sparse (50k nnz): $(format_bytes(alloc2))")

println()

println("3.2 NTuple construction methods")
# Test the difference between generator and ntuple
nΛder, nXder, nUder = 2, 3, 1
state_Λ = (zeros(100), zeros(100))
state_X = (zeros(100), zeros(100), zeros(100))
state_U = (zeros(50),)
index_X = [1, 2, 3, 4, 5]

# Current method (generator)
function tuple_gen(state_Λ, index_X, nΛder)
    NTuple{nΛder}(λ[index_X] for λ∈state_Λ)
end

# Proposed method (ntuple)
function tuple_ntuple(state_Λ, index_X, nΛder)
    ntuple(i -> state_Λ[i][index_X], Val(nΛder))
end

# Warm up
tuple_gen(state_Λ, index_X, nΛder)
tuple_ntuple(state_Λ, index_X, nΛder)

alloc_gen = @allocated tuple_gen(state_Λ, index_X, nΛder)
alloc_ntuple = @allocated tuple_ntuple(state_Λ, index_X, nΛder)

println("  Generator expression alloc: $(format_bytes(alloc_gen))")
println("  ntuple with Val alloc: $(format_bytes(alloc_ntuple))")

println()

# ============================================================================
# Section 4: Element-level benchmarks (requires full model setup)
# ============================================================================
println("=" ^ 70)
println("Section 4: Element assembly benchmarks")
println("=" ^ 70)
println()
println("  (Skipped - requires full model setup. Use existing tests for validation.)")
println()

# ============================================================================
# Section 5: Full solver benchmark (if time permits)
# ============================================================================
println("=" ^ 70)
println("Section 5: Summary")
println("=" ^ 70)
println()

println("Benchmark Results Summary:")
println("-" ^ 70)
for r in results
    @printf "  %-35s: %8.3f ms (allocs: %d)\n" r.name r.mean r.allocs
end

println()
println("=" ^ 70)
println("Benchmark complete")
println("=" ^ 70)

# Return results for potential comparison
results