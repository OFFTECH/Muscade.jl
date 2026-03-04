using Pkg
Pkg.activate(dirname(@__DIR__))
Pkg.instantiate()

using Muscade
using SparseArrays
using BenchmarkTools
using Printf

println("=" ^ 70)
println("PERFORMANCE COMPARISON: Before vs After SIMD Optimization")
println("=" ^ 70)
println()

# OLD implementation (without SIMD)
function zero_old!(out::AbstractArray)
    for i∈eachindex(out)
        out[i] = 0
    end
end
function zero_old!(out::AbstractSparseArray)
    for i∈eachindex(out.nzval)
        out.nzval[i] = 0
    end
end

# NEW implementation (with SIMD - from current Muscade)
# Uses Muscade.zero! which has @simd and @inbounds

println("=" ^ 70)
println("TEST 1: Dense Vector zero! (n = 10,000)")
println("=" ^ 70)

v_old = zeros(10_000)
v_new = zeros(10_000)

b_old = @benchmark zero_old!($v_old) samples=1000 evals=10
b_new = @benchmark Muscade.zero!($v_new) samples=1000 evals=10

println("\nOLD (no SIMD):")
display(b_old)
println("\n\nNEW (with @simd @inbounds):")
display(b_new)

speedup_1 = median(b_old.times) / median(b_new.times)
println("\n\nSpeedup: $(round(speedup_1, digits=2))x")

println()
println("=" ^ 70)
println("TEST 2: Dense Vector zero! (n = 100,000)")
println("=" ^ 70)

v_old2 = zeros(100_000)
v_new2 = zeros(100_000)

b_old2 = @benchmark zero_old!($v_old2) samples=1000 evals=10
b_new2 = @benchmark Muscade.zero!($v_new2) samples=1000 evals=10

println("\nOLD (no SIMD):")
display(b_old2)
println("\n\nNEW (with @simd @inbounds):")
display(b_new2)

speedup_2 = median(b_old2.times) / median(b_new2.times)
println("\n\nSpeedup: $(round(speedup_2, digits=2))x")

println()
println("=" ^ 70)
println("TEST 3: Sparse Matrix zero! (1000x1000, ~5000 nnz)")
println("=" ^ 70)

sp_old = sprand(1000, 1000, 0.005)
sp_new = copy(sp_old)

b_old3 = @benchmark zero_old!($sp_old) samples=1000 evals=10
b_new3 = @benchmark Muscade.zero!($sp_new) samples=1000 evals=10

println("\nOLD (no SIMD):")
display(b_old3)
println("\n\nNEW (with @simd @inbounds):")
display(b_new3)

speedup_3 = median(b_old3.times) / median(b_new3.times)
println("\n\nSpeedup: $(round(speedup_3, digits=2))x")

println()
println("=" ^ 70)
println("TEST 4: Sparse Matrix zero! (10000x10000, ~50000 nnz)")
println("=" ^ 70)

sp_old4 = sprand(10000, 10000, 0.0005)
sp_new4 = copy(sp_old4)

b_old4 = @benchmark zero_old!($sp_old4) samples=1000 evals=10
b_new4 = @benchmark Muscade.zero!($sp_new4) samples=1000 evals=10

println("\nOLD (no SIMD):")
display(b_old4)
println("\n\nNEW (with @simd @inbounds):")
display(b_new4)

speedup_4 = median(b_old4.times) / median(b_new4.times)
println("\n\nSpeedup: $(round(speedup_4, digits=2))x")

println()
println("=" ^ 70)
println("SUMMARY")
println("=" ^ 70)
println()
println("| Test                          | Old (μs) | New (μs) | Speedup |")
println("|-------------------------------|----------|----------|---------|")
@printf "| Dense 10k                     | %8.1f | %8.1f | %5.2fx  |\n" median(b_old.times)/1000 median(b_new.times)/1000 speedup_1
@printf "| Dense 100k                    | %8.1f | %8.1f | %5.2fx  |\n" median(b_old2.times)/1000 median(b_new2.times)/1000 speedup_2
@printf "| Sparse 5k nnz                 | %8.2f | %8.2f | %5.2fx  |\n" median(b_old3.times)/1000 median(b_new3.times)/1000 speedup_3
@printf "| Sparse 50k nnz                | %8.1f | %8.1f | %5.2fx  |\n" median(b_old4.times)/1000 median(b_new4.times)/1000 speedup_4

println()
println("Note: Times in microseconds (μs)")
println("=" ^ 70)
