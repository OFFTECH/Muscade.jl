"""
Performance Regression Tests for Muscade.jl

This script runs benchmarks and compares against baseline performance metrics.
If performance degrades beyond a threshold, the test fails.

Usage:
    julia --project=. benchmark/bench_regression.jl [--baseline baseline.jld2] [--threshold 1.2]

Arguments:
    --baseline FILE     File containing baseline performance metrics (JLD2 format)
    --threshold FLOAT    Maximum allowed performance degradation (default: 1.2 = 20% slowdown allowed)
    --update FILE      Update baseline file with current performance
"""

using Pkg
Pkg.activate(".")
using Muscade
using BenchmarkTools

# Include test elements
include("../test/SomeElements.jl")

# Simple file-based storage without external dependencies
"""
Save benchmark results to a simple text file.
"""
function save_results_simple(results, filename)
    open(filename, "w") do f
        for r in results
            println(f, "$(r["name"]),$(r["nelements"]),$(r["nsteps"]),",
                    "$(r["min_time_ms"]),$(r["median_time_ms"]),",
                    "$(r["memory_mb"]),$(r["threads"])")
        end
    end
    println("Results saved to $filename")
end

"""
Load benchmark results from a simple text file.
"""
function load_results_simple(filename)
    if !isfile(filename)
        println("Baseline file not found: $filename")
        return []
    end

    results = []
    open(filename, "r") do f
        for line in eachline(f)
            if isempty(strip(line))
                continue
            end
            parts = split(strip(line), ',')
            if length(parts) >= 7
                push!(results, Dict{String, Any}(
                    "name" => parts[1],
                    "nelements" => parse(Int, parts[2]),
                    "nsteps" => parse(Int, parts[3]),
                    "min_time_ms" => parse(Float64, parts[4]),
                    "median_time_ms" => parse(Float64, parts[5]),
                    "memory_mb" => parse(Float64, parts[6]),
                    "threads" => parse(Int, parts[7])
                ))
            end
        end
    end
    return results
end

"""
Create a benchmark model with SdofOscillator elements.

Args:
    nelements: Number of SdofOscillator elements to create

Returns:
    model: The benchmark model
    state0: Initial state
"""
function create_benchmark_model(nelements)
    model = Model(:BenchModel)
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
Run a single benchmark case.

Args:
    name: Name of the benchmark case
    nelements: Number of elements
    nsteps: Number of time steps

Returns:
    Dictionary with benchmark results
"""
function run_benchmark_case(name, nelements, nsteps)
    println("Running benchmark: $name (elements=$nelements, steps=$nsteps)")

    model, state0 = create_benchmark_model(nelements)

    # Warm-up run
    try
        solve(SweepX{2}; initialstate=deepcopy(state0), time=0.01:0.1:(0.01 * nsteps), verbose=false)
    catch
        # Warm-up might fail for some configurations, that's OK
    end

    # Benchmark run
    b = @benchmark begin
        solve(SweepX{2}; initialstate=deepcopy($state0),
              time=0.01:0.1:(0.01 * $nsteps), verbose=false)
    end samples=10 evals=1

    # Extract key metrics
    min_time = minimum(b).time / 1e6  # Convert ns to ms
    median_time = median(b).time / 1e6
    memory_mb = BenchmarkTools.memory(b) / 1024 / 1024

    result = Dict{String, Any}(
        "name" => name,
        "nelements" => nelements,
        "nsteps" => nsteps,
        "min_time_ms" => min_time,
        "median_time_ms" => median_time,
        "memory_mb" => memory_mb,
        "threads" => Threads.nthreads()
    )

    println("  Median: $(round(median_time, digits=2)) ms")
    println("  Min:    $(round(min_time, digits=2)) ms")
    println("  Memory: $(round(memory_mb, digits=2)) MiB")

    return result
end

"""
Run all benchmark cases.

Returns:
    Array of benchmark result dictionaries
"""
function run_all_benchmarks()
    println("=" ^ 60)
    println("Performance Regression Benchmarks")
    println("=" ^ 60)
    println("Julia version: $(VERSION)")
    println("Threads: $(Threads.nthreads())")
    println()

    results = []

    # Small model - quick test
    push!(results, run_benchmark_case("small_100_elements", 100, 20))

    # Medium model
    push!(results, run_benchmark_case("medium_500_elements", 500, 20))

    # Large model - tests threading
    push!(results, run_benchmark_case("large_1000_elements", 1000, 20))

    return results
end

"""
Compare current results against baseline.

Args:
    results: Current benchmark results
    baseline: Baseline benchmark results
    threshold: Maximum allowed degradation ratio

Returns:
    true if all metrics are within threshold, false otherwise
"""
function compare_against_baseline(results, baseline, threshold)
    println()
    println("=" ^ 60)
    println("Comparison against baseline")
    println("=" ^ 60)

    # Create lookup by name for baseline
    baseline_dict = Dict(b["name"] => b for b in baseline)
    all_passed = true

    for result in results
        name = result["name"]
        if !haskey(baseline_dict, name)
            println("  $name: NO BASELINE (new benchmark)")
            continue
        end

        b = baseline_dict[name]
        current_median = result["median_time_ms"]
        baseline_median = b["median_time_ms"]

        ratio = current_median / baseline_median
        diff_pct = (ratio - 1) * 100

        # Check if configuration changed (different threads)
        if result["threads"] != b["threads"]
            println("  $name: SKIPPED (thread count changed: $(b["threads"]) -> $(result["threads"]))")
            continue
        end

        # Check threshold
        if ratio > threshold
            println("  $name: FAILED")
            println("    Baseline: $(round(baseline_median, digits=2)) ms")
            println("    Current:  $(round(current_median, digits=2)) ms")
            println("    Degradation: $(round(diff_pct, digits=1))% (threshold: $(round((threshold - 1) * 100, digits=0))%)")
            all_passed = false
        else
            println("  $name: PASSED ($(round(diff_pct, digits=1))% change)")
        end
    end

    return all_passed
end

"""
Main function to run performance regression tests.

Returns:
    exit code (0 for success, 1 for failure)
"""
function main()
    # Parse command line arguments manually
    baseline_file = "benchmark/baseline.txt"
    threshold = 1.2
    update_baseline = false

    for i in 1:length(ARGS)
        if ARGS[i] == "--baseline" && i < length(ARGS)
            baseline_file = ARGS[i+1]
        elseif ARGS[i] == "--threshold" && i < length(ARGS)
            threshold = parse(Float64, ARGS[i+1])
        elseif ARGS[i] == "--update" && i < length(ARGS)
            baseline_file = ARGS[i+1]
            update_baseline = true
        elseif ARGS[i] == "-b" && i < length(ARGS)
            baseline_file = ARGS[i+1]
        elseif ARGS[i] == "-t" && i < length(ARGS)
            threshold = parse(Float64, ARGS[i+1])
        elseif ARGS[i] == "-u" && i < length(ARGS)
            baseline_file = ARGS[i+1]
            update_baseline = true
        end
    end

    # Run benchmarks
    results = run_all_benchmarks()

    # Handle --update flag
    if update_baseline
        save_results_simple(results, baseline_file)
        return 0
    end

    # Load and compare against baseline
    baseline = load_results_simple(baseline_file)

    if isempty(baseline)
        println()
        println("No baseline found. Creating baseline file...")
        save_results_simple(results, baseline_file)
        println("Run with --baseline $baseline_file to compare in future runs.")
        return 0
    end

    # Compare and report
    passed = compare_against_baseline(results, baseline, threshold)

    println()
    if passed
        println("All performance regression tests PASSED")
        return 0
    else
        println("Some performance regression tests FAILED")
        println()
        println("If the degradation is expected (e.g., new features), update baseline with:")
        println("  julia --project=. benchmark/bench_regression.jl --update $baseline_file")
        return 1
    end
end

# Run main function and exit with appropriate code
exit(main())
