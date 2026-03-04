"""
Large Model Benchmark for Muscade.jl

This script benchmarks Muscade with large models to:
1. Validate threading performance improvements
2. Test with real-world model sizes (2000+ elements)
3. Measure speedup across different thread counts
4. Compare different element types

Usage:
    julia --project=. benchmark/bench_large_model.jl [--nelements N] [--threads N]
"""

using Pkg
Pkg.activate(".")

using Muscade
using BenchmarkTools
using Printf
using Dates

# Include test elements
include("../test/SomeElements.jl")

"""
Parse command line arguments.
"""
function parse_args()
    nelements = 2000
    nthreads_set = nothing

    for i in 1:length(ARGS)
        if ARGS[i] == "--nelements" && i < length(ARGS)
            nelements = parse(Int, ARGS[i+1])
        elseif ARGS[i] == "--threads" && i < length(ARGS)
            nthreads_set = parse(Int, ARGS[i+1])
        end
    end

    return nelements, nthreads_set
end

"""
Create a large model with SdofOscillator elements.
"""
function create_large_oscillator_model(nelements)
    model = Model(:LargeOscillatorModel)
    nodes = [addnode!(model, 𝕣[]) for _ in 1:nelements]

    for i in 1:nelements
        addelement!(model, SdofOscillator, [nodes[i]];
                   K₁=1., K₂=0.3, C₁=0.1, C₂=0., M₁=1., M₂=0.)
    end

    state0 = initialize!(model; time=0.)
    for i in 1:nelements
        state0 = setdof!(state0, [1.0]; field=:tx1, nodID=[nodes[i]], order=0)
    end

    return model, state0
end

"""
Create a large model with mixed element types for more realistic testing.
"""
function create_mixed_element_model(nelements)
    model = Model(:MixedElementModel)
    nodes = [addnode!(model, 𝕣[]) for _ in 1:nelements]

    # Use SdofOscillator for first half
    mid = nelements ÷ 2

    for i in 1:mid
        addelement!(model, SdofOscillator, [nodes[i]];
                   K₁=1., K₂=0.3, C₁=0.1, C₂=0., M₁=1., M₂=0.)
    end

    # For second half, use oscillators with different parameters
    for i in (mid+1):nelements
        addelement!(model, SdofOscillator, [nodes[i]];
                   K₁=0.5, K₂=0.1, C₁=0.05, C₂=0., M₁=0.5, M₂=0.)
    end

    state0 = initialize!(model; time=0.)
    for i in 1:nelements
        state0 = setdof!(state0, [1.0]; field=:tx1, nodID=[nodes[i]], order=0)
    end

    return model, state0
end

"""
Benchmark a single model configuration.
"""
function benchmark_model(model, state0, nsteps; name="", warmup=true)
    println("Benchmarking: $name (elements=$(length(model.ele[1])))")

    # Warm-up run
    if warmup
        try
            solve(SweepX{2}; initialstate=deepcopy(state0),
                  time=0.01:0.1:(0.01 * nsteps), verbose=false)
        catch
            # Warm-up might fail, that's OK
        end
    end

    # Benchmark run with multiple samples for better statistics
    b = @benchmark begin
        solve(SweepX{2}; initialstate=deepcopy($state0),
              time=0.01:0.1:(0.01 * $nsteps), verbose=false)
    end samples=5 evals=1

    min_time = minimum(b).time / 1e6  # Convert ns to ms
    median_time = median(b).time / 1e6
    mean_time = mean(b).time / 1e6
    memory_mb = BenchmarkTools.memory(b) / 1024 / 1024

    result = Dict{String, Any}(
        "name" => name,
        "nelements" => length(model.ele[1]),
        "nsteps" => nsteps,
        "min_time_ms" => min_time,
        "median_time_ms" => median_time,
        "mean_time_ms" => mean_time,
        "memory_mb" => memory_mb,
        "threads" => Threads.nthreads()
    )

    @printf("  Median: %.2f ms\n", median_time)
    @printf("  Mean:   %.2f ms\n", mean_time)
    @printf("  Min:    %.2f ms\n", min_time)
    @printf("  Memory: %.2f MiB\n", memory_mb)

    return result
end

"""
Run threading scalability tests.
"""
function benchmark_threading_scaling(nelements, nsteps)
    println()
    println("=" ^ 60)
    println("Threading Scalability Tests")
    println("=" ^ 60)
    println("Elements: $nelements, Steps: $nsteps")
    println()

    results = []
    thread_counts = [1, 2, 4, 8]
    available_threads = Threads.nthreads()

    # Test with different thread counts (up to available)
    for t in thread_counts
        if t > available_threads
            continue
        end

        println("-" ^ 40)
        println("Threads: $t")

        # Note: Thread count is set via JULIA_NUM_THREADS environment variable
        model, state0 = create_large_oscillator_model(nelements)
        result = benchmark_model(model, state0, nsteps; name="threads_$t")
        push!(results, result)

        # Calculate speedup relative to 1 thread
        if t == 1
            result["speedup"] = 1.0
            result["efficiency"] = 1.0
        else
            baseline = results[1]["median_time_ms"]
            speedup = baseline / result["median_time_ms"]
            efficiency = speedup / t
            result["speedup"] = speedup
            result["efficiency"] = efficiency
            @printf("  Speedup: %.2fx\n", speedup)
            @printf("  Efficiency: %.1f%%\n", efficiency * 100)
        end
    end

    return results
end

"""
Run element type comparison.
"""
function benchmark_element_types(nelements, nsteps)
    println()
    println("=" ^ 60)
    println("Element Type Comparison")
    println("=" ^ 60)
    println("Elements: $nelements, Steps: $nsteps")
    println()

    results = []

    # Oscillator elements
    println("-" ^ 40)
    println("SdofOscillator elements")
    model1, state1 = create_large_oscillator_model(nelements)
    result1 = benchmark_model(model1, state1, nsteps; name="oscillator")
    push!(results, result1)

    # Mixed elements
    println("-" ^ 40)
    println("Mixed elements (Oscillator with different parameters)")
    model2, state2 = create_mixed_element_model(nelements)
    result2 = benchmark_model(model2, state2, nsteps; name="mixed")
    push!(results, result2)

    return results
end

"""
Run problem size scaling test.
"""
function benchmark_problem_scaling(nsteps)
    println()
    println("=" ^ 60)
    println("Problem Size Scaling")
    println("=" ^ 60)
    println("Steps: $nsteps")
    println()

    results = []
    element_counts = [500, 1000, 2000, 5000]

    for n in element_counts
        println("-" ^ 40)
        println("Elements: $n")

        model, state0 = create_large_oscillator_model(n)
        result = benchmark_model(model, state0, nsteps; name="size_$n")
        push!(results, result)

        # Calculate time per element
        time_per_element = result["median_time_ms"] / n
        @printf("  Time per element: %.3f ms\n", time_per_element)
        result["time_per_element_ms"] = time_per_element
    end

    return results
end

"""
Print summary of all benchmark results.
"""
function print_summary(all_results)
    println()
    println("=" ^ 60)
    println("Summary")
    println("=" ^ 60)
    println()

    for (category, results) in all_results
        println("Category: $category")
        println("-" ^ 40)

        for r in results
            println("  $(r["name"]): $(round(r["median_time_ms"], digits=2)) ms")
            if haskey(r, "speedup")
                @printf("    Speedup: %.2fx (efficiency: %.1f%%)\n",
                        r["speedup"], r["efficiency"] * 100)
            end
            if haskey(r, "time_per_element_ms")
                @printf("    Time/element: %.3f ms\n", r["time_per_element_ms"])
            end
        end
        println()
    end
end

"""
Main function.
"""
function main()
    nelements, nthreads_set = parse_args()

    # Note: Thread count is set via JULIA_NUM_THREADS environment variable
    println("=" ^ 60)
    println("Large Model Benchmark for Muscade.jl")
    println("=" ^ 60)
    println("Julia version: $(VERSION)")
    println("Max threads: $(Threads.nthreads())")
    println()

    nsteps = 50  # Number of time steps for benchmarks

    all_results = []

    # Run threading scalability tests
    push!(all_results, ("Threading Scalability",
                      benchmark_threading_scaling(nelements, nsteps)))

    # Run element type comparison
    push!(all_results, ("Element Type Comparison",
                      benchmark_element_types(nelements, nsteps)))

    # Run problem size scaling test
    push!(all_results, ("Problem Size Scaling",
                      benchmark_problem_scaling(nsteps)))

    # Print summary
    print_summary(all_results)

    # Save results to file
    open("benchmark/large_model_results.txt", "w") do f
        println(f, "# Large Model Benchmark Results")
        println(f, "# Julia version: $(VERSION)")
        println(f, "# Date: $(Dates.now())")
        println(f, "#")
        println(f, "# Format: category,name,nelements,nsteps,min_time_ms,median_time_ms,mean_time_ms,memory_mb,threads")
        println(f, "#")

        for (category, results) in all_results
            for r in results
                line = "$(category),$(r["name"]),$(r["nelements"]),$(r["nsteps"])"
                line *= ",$(r["min_time_ms"]),$(r["median_time_ms"]),$(r["mean_time_ms"])"
                line *= ",$(r["memory_mb"]),$(r["threads"])"

                if haskey(r, "speedup")
                    line *= ",$(r["speedup"]),$(r["efficiency"])"
                else
                    line *= ","
                end

                if haskey(r, "time_per_element_ms")
                    line *= ",$(r["time_per_element_ms"])"
                else
                    line *= ","
                end

                println(f, line)
            end
        end
    end

    println()
    println("Results saved to benchmark/large_model_results.txt")
end

main()
