using Pkg
Pkg.activate(dirname(@__DIR__))
Pkg.instantiate()

using Muscade
using SparseArrays
using BenchmarkTools
using Printf

println("=" ^ 60)
println("Benchmark: zero! function performance")
println("=" ^ 60)
println()

# Test 1: Dense array zero!
println("1. Dense Array zero! (Vector{Float64}, n=10_000)")
v = zeros(10_000)
b1 = @benchmark Muscade.zero!($v) samples=1000 evals=10
display(b1)
println()

# Test 2: Large dense array
println("2. Dense Array zero! (Vector{Float64}, n=100_000)")
v2 = zeros(100_000)
b2 = @benchmark Muscade.zero!($v2) samples=1000 evals=10
display(b2)
println()

# Test 3: Sparse array zero!
println("3. Sparse Array zero! (1000x1000, ~5000 non-zeros)")
sp = sprand(1000, 1000, 0.005)
b3 = @benchmark Muscade.zero!($sp) samples=1000 evals=10
display(b3)
println()

# Test 4: Large sparse array
println("4. Sparse Array zero! (10000x10000, ~50000 non-zeros)")
sp2 = sprand(10000, 10000, 0.0005)
b4 = @benchmark Muscade.zero!($sp2) samples=1000 evals=10
display(b4)
println()

# Test 5: Matrix
println("5. Dense Matrix zero! (100x100)")
m = zeros(100, 100)
b5 = @benchmark Muscade.zero!($m) samples=1000 evals=10
display(b5)
println()

println()
println("=" ^ 60)
println("Memory allocation test")
println("=" ^ 60)
println()

# Memory allocation test
v_alloc = zeros(100_000)
alloc = @allocated Muscade.zero!(v_alloc)
println("Allocations for zero! on 100k vector: $alloc bytes")

sp_alloc = sprand(10000, 10000, 0.0005)
alloc_sp = @allocated Muscade.zero!(sp_alloc)
println("Allocations for zero! on sparse (50k nnz): $alloc_sp bytes")

println()
println("=" ^ 60)
println("Benchmark complete")
println("=" ^ 60)
