# Cache Efficiency Profiling Guide for Muscade.jl

using Pkg
Pkg.activate(".")

using Muscade
using Profile
using BenchmarkTools
using Printf

println("=" ^60)
println("CACHE EFFICIENCY PROFILING GUIDE")
println("=" ^ 60)

# ============================================================
# Method 1: Built-in Julia Profiler
# ============================================================

"""
Basic profiling using Julia's built-in profiler.
Shows where time is spent, which helps identify cache issues.
"""
function profile_basic()
    println("\n--- Method 1: Built-in Profiler ---")
    
    include("../test/SomeElements.jl")
    
    # Create model
    model = Model(:ProfileModel)
    nodes = [addnode!(model, 𝕣[]) for _ in 1:500]
    for i in 1:500
        addelement!(model, SdofOscillator, [nodes[i]]; K₁=1., K₂=0.3, C₁=0.1, C₂=0., M₁=1., M₂=0.)
    end
    
    state = initialize!(model; time=0.)
    for i in 1:500
        state = setdof!(state, [1.0]; field=:tx1, nodID=[nodes[i]], order=0)
    end
    
    # Warmup
    solve(SweepX{2}; initialstate=deepcopy(state), time=0.01:0.1:0.11, verbose=false)
    
    # Profile
    Profile.clear()
    @profile for _ in 1:10
        solve(SweepX{2}; initialstate=deepcopy(state), time=0.01:0.1:0.2, verbose=false)
    end
    
    # Print results
    println("\nTop functions by time:")
    Profile.print(format=:flat, sortedby=:count, maxdepth=15)
    
    return model, state
end

# ============================================================
# Method 2: Memory Allocation Profiling
# ============================================================

"""
Profile memory allocations to identify cache pressure from GC.
"""
function profile_allocations(model, state)
    println("\n--- Method 2: Allocation Profiling ---")
    
    # Use @allocated to measure per-operation allocations
    allocations = @allocated solve(SweepX{2}; initialstate=deepcopy(state), time=0.01:0.1:0.2, verbose=false)
    println("Total allocations per solve: $(allocations / 1024 / 1024) MB")
    
    # Detailed allocation tracking
    println("\nAllocation breakdown:")
    println("  deepcopy(state): $(@allocated(deepcopy(state)) / 1024) KB")
end

# ============================================================
# Method 3: Time Breakdown by Function
# ============================================================

"""
Measure time spent in each major component.
"""
function profile_components(model, state)
    println("\n--- Method 3: Component Timing ---")
    
    # Get assembler
    dis = state.dis
    out, asm, Xdofgr = Muscade.prepare(AssemblySweepX{2}, model, dis)
    
    # Profile assembly
    b_assembly = @benchmark assemble!{:iter}($out, $asm, $dis, $model, $state, 0.01, (;)) samples=100
    println("Assembly time: $(round(median(b_assembly).time / 1e6, digits=3)) ms")
    
    # Profile LU factorization
    assemble!{:iter}(out, asm, dis, model, state, 0.01, (;))
    b_lu = @benchmark lu($out.Lλx) samples=100
    println("LU factorization: $(round(median(b_lu).time / 1e6, digits=3)) ms")
    
    # Profile solve
    F = lu(out.Lλx)
    b_solve = @benchmark ldiv!($(similar(out.Lλ)), $F, $out.Lλ) samples=100
    println("Linear solve: $(round(median(b_solve).time / 1e6, digits=3)) ms")
    
    # Calculate percentages
    total = median(b_assembly).time + median(b_lu).time + median(b_solve).time
    println("\nTime breakdown:")
    @printf("  Assembly:        %5.1f%%\n", 100 * median(b_assembly).time / total)
    @printf("  LU factorization: %5.1f%%\n", 100 * median(b_lu).time / total)
    @printf("  Linear solve:    %5.1f%%\n", 100 * median(b_solve).time / total)
end

# ============================================================
# Method 4: Cache Miss Analysis (Linux only)
# ============================================================

"""
Instructions for using perf on Linux to measure cache misses.
"""
function cache_miss_instructions()
    println("\n--- Method 4: Cache Miss Analysis (Linux) ---")
    println("""
    On Linux, use perf to measure cache misses:
    
    # Run Julia with perf:
    perf stat -e cache-references,cache-misses,L1-dcache-loads,L1-dcache-load-misses \\
        julia --threads=4 benchmark/bench_hpc_final.jl
    
    # For detailed analysis:
    perf record -g julia --threads=4 benchmark/bench_hpc_final.jl
    perf report
    
    # Key metrics to watch:
    # - cache-misses / cache-references (should be < 10%)
    # - L1-dcache-load-misses / L1-dcache-loads (should be < 5%)
    """)
end

# ============================================================
# Method 5: Memory Access Pattern Analysis
# ============================================================

"""
Analyze memory access patterns in assembly.
"""
function analyze_access_patterns()
    println("\n--- Method 5: Memory Access Pattern Analysis ---")
    
    println("""
    Key questions for cache efficiency:
    
    1. ELEMENT ACCESS PATTERN:
       - Elements are accessed sequentially: GOOD for cache
       - Each element accesses its own nodes: may have poor locality
       
    2. NODE DATA ACCESS:
       - Node coordinates accessed per element
       - Random access if elements not spatially ordered
       - Consider: element reordering by spatial locality
    
    3. ASSEMBLY MATRIX ACCESS:
       - Sparse matrix assembly via nzval
       - Indirect indexing: asm[k, iele] where k varies
       - Non-contiguous writes to global arrays
    
    4. STATE VECTOR ACCESS:
       - state.X[ider][index.X] - indirect indexing
       - May cause cache misses if index.X is not contiguous
    
    TO PROFILE:
    
    a) Add timing around specific loops:
       @timeit to "element loop" for iele = 1:nele
           @timeit to "index lookup" index = dis.index[iele]
           @timeit to "residual" R = getresidual(...)
           @timeit to "assembly" add_∂!(...)
       end
    
    b) Use TimerOutputs.jl for nested timing:
       using TimerOutputs
       to = TimerOutput()
       # ... code with @timeit ...
       display(to)
    """)
end

# ============================================================
# Method 6: Recommended Profiling Approach
# ============================================================

function recommended_approach()
    println("\n" * "=" ^60)
    println("RECOMMENDED PROFILING APPROACH")
    println("=" ^ 60)
    
    println("""
    STEP 1: Identify Hot Spots
    --------------------------
    julia> using Profile
    julia> @profile solve(...)  # Run 10-100 times
    julia> Profile.print(format=:flat, sortedby=:count)
    
    Look for functions with high % time.
    
    STEP 2: Measure Component Times
    -------------------------------
    Use BenchmarkTools to time individual components:
    - Assembly time
    - LU factorization time
    - Linear solve time
    
    If assembly > 50%, focus there.
    
    STEP 3: Analyze Memory Access (Linux)
    -------------------------------------
    perf stat -e cache-misses,cache-references \\
        julia --threads=4 your_benchmark.jl
    
    If cache-miss rate > 10%, consider reordering.
    
    STEP 4: Add Instrumentation
    ---------------------------
    using TimerOutputs
    
    const TO = TimerOutput()
    
    function assemble_!(...)
        @timeit TO "element loop" begin
            for iele = 1:nele
                @timeit TO "get index" index = dis.index[iele]
                @timeit TO "residual" R = getresidual(...)
                @timeit TO "assemble" add_∂!(...)
            end
        end
    end
    
    # After solve:
    display(TO)
    
    STEP 5: Compare Strategies
    --------------------------
    - Baseline timing
    - After element reordering
    - After buffer reorganization
    - After prefetching
    
    DOCUMENT RESULTS in HPC_IMPROVEMENTS.md
    """)
end

# ============================================================
# Run All
# ============================================================

println("\nRun these functions to profile:")
println("  profile_basic()           - Basic time profiling")
println("  profile_allocations()     - Memory allocation analysis")
println("  profile_components()      - Component timing breakdown")
println("  cache_miss_instructions() - Linux perf instructions")
println("  analyze_access_patterns() - Access pattern analysis")
println("  recommended_approach()    - Full profiling guide")

recommended_approach()
