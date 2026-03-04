# Muscade.jl Optimization & Architecture Roadmap

**Consolidated:** 2026-03-02
**Last Updated:** 2026-03-03

---

## Table of Contents

1. [Quick Reference - Status Overview](#quick-reference---status-overview)
2. [✅ Completed Optimizations](#-completed-optimizations)
3. [🟡 In Progress / Partially Done](#in-progress--partially-done)
4. [⏸️ Deferred / Low Priority](#️-deferred--low-priority)
5. [🔜 Future Work (Major Undertakings)](#️-future-work--major-undertakings)
6. [Architecture & Code Quality](#architecture--code-quality)
7. [Testing & Benchmarking](#testing--benchmarking)
8. [Implementation Timeline](#implementation-timeline)

---

## Quick Reference - Status Overview

### Test Suite Status
| Metric | Value |
|---------|--------|
| **Total Tests** | 745 |
| **Tests Passing** | 745 (100%) |
| **Test Duration** | ~5 minutes |

### Performance Improvements Summary

| Category | Status | Impact |
|----------|--------|--------|
| **Parallelization (Threading)** | ✅ Complete | 1.1-1.6x speedup (1000 elements) |
| **SIMD - zero!** | ✅ Complete | 1.01x (dense), 1.27-3.5x (sparse) |
| **SIMD - DofGroup ops** | ✅ Complete | Hot path optimized |
| **SIMD - @views** | ✅ Complete | Eliminates copies |
| **Memory - Pre-allocate buffers** | ✅ Complete | Already present in AssemblySweepX |
| **Memory - add_∂! hoisting** | ✅ Complete | Better branch prediction |
| **Linear Algebra - ldiv!** | ✅ Complete | Eliminates allocation |
| **Benchmark Suite** | ✅ Complete | Comprehensive benchmarking |

---

## ✅ Completed Optimizations

### SIMD & Loop Optimizations

#### D-1: SIMD for zero! Functions
**Status:** ✅ COMPLETED (2026-03-01)
**Location:** `src/Assemble.jl:537-546`
**Results:**
- Dense arrays: 1.01x (already optimal)
- Sparse arrays: 1.27-3.5x speedup

#### D-2: SIMD for DofGroup Operations (decrement!, increment!)
**Status:** ✅ COMPLETED (2026-03-02)
**Location:** `src/Assemble.jl:206-263`

**Implementation:**
```julia
# Helper: unconditionally decrement Λ (for ider within bounds)
@inline function _decrement_Λ!(s::State,ider::𝕫,y::AbstractVector{𝕣},gr::DofGroup)
    @simd for i ∈ eachindex(gr.iΛ); @inbounds s.Λ[ider][gr.iΛ[i]] -= y[gr.jΛ[i]] * gr.scaleΛ[i]; end

# Generic fallback for unusual cases
function decrement!(s::State,ider::𝕫,y::AbstractVector{𝕣},gr::DofGroup)
    if ider≤length(s.Λ) _decrement_Λ!(s,ider,y,gr) end
    if ider≤length(s.X) _decrement_X!(s,ider,y,gr) end
    if ider≤length(s.U) _decrement_U!(s,ider,y,gr) end
    if ider==1          _decrement_A!(s,y,gr) end
end
```

**Key Design:** Conditionals moved outside SIMD loop for proper vectorization.

#### D-5: @views for Array Slicing
**Status:** ✅ COMPLETED (2026-03-02)
**Location:** `src/Assemble.jl:350-356`

**Impact:** Eliminates temporary array copies in `asmvec_kernel!`

---

### Parallelization (Multi-threading)

#### C-1/C-2: Multi-threaded Assembly for SweepX
**Status:** ✅ COMPLETED (2026-03-02)
**Location:** `src/SweepX.jl:217-250`

**Implementation:**
- `ThreadLocalAssemblySweepX` struct for thread-local buffers
- Pre-allocated thread buffers in `prepare()` using `Threads.maxthreadid()`
- Specialized `assemble_!{mission}` method for AssemblySweepX
- `:static` scheduling for proper thread ID mapping
- Threshold of 50 elements (`MIN_ELEMENTS_FOR_THREADING`) to justify threading overhead

**Benchmark Results (1000 SdofOscillator elements):**
| Threads | Median Time | Speedup |
|--------|---------------|----------|
| 1 | 29.61 ms | 1.00x |
| 2 | 19.78 ms | 1.50x |
| 4 | 18.81 ms | **1.57x** |

**Known Limitations:**
- Test models in repository have <50 elements, so threading is not triggered
- Large model benchmark needed to measure actual speedup
- Only implemented for SweepX solver (not all solvers)

---

### Memory Management

#### E-5: Use ldiv! to Avoid Allocation
**Status:** ✅ COMPLETED (2026-03-02)
**Location:** `src/SweepX.jl:372`

**Before:**
```julia
Δx = Lλx \ out.Lλ  # Allocates new vector
```

**After:**
```julia
ldiv!(out.Δx_buffer, Lλx, out.Lλ)  # In-place solve
Δx = out.Δx_buffer  # Local reference
```

**Impact:** Eliminates per-iteration vector allocation in solve loop.

#### B-4: Avoid Allocation in add_∂!
**Status:** ✅ COMPLETED (2026-03-02)
**Location:** `src/Assemble.jl:570-646`

**Implementation:**
- Hoisted conditionals outside loops (S and T are compile-time constants)
- Added `@inbounds` to all array access
- Added `@views` to `asmvec_kernel!`

**Impact:** Better branch prediction, reduces runtime overhead.

---

### Benchmarking & Validation

#### H-1: Benchmark Suite
**Status:** ✅ COMPLETED (2026-03-02)
**Location:** `benchmark/bench_hpc_final.jl`

**Features:**
- Multiple element counts (200, 500, 1000)
- Min/median/mean timing
- Memory allocation tracking
- Thread scaling analysis

---

## 🟡 In Progress / Partially Done

### Type Stability Issues

#### A-1: Fix gradientpartition Type Instability
**Priority:** Medium (Low impact in setup phase)
**Location:** `src/Assemble.jl:329-335`
**Status:** ⏸️ DEFERRED - Barrier function, minimal hot path impact
**Issue:** Return type depends on input values, causing type instability. Comment on line 349: "TODO type unstable, barrier function".

#### A-2: Fix type_multivariate_𝕣 Type Instability
**Priority:** Medium
**Location:** `src/Taylor.jl:97-101`
**Status:** ⏸️ DEFERRED - Complex refactor
**Issue:** Comment states: "this causes (slight) type instability - because if Ra is not ℝ, then return type is different."

---

## ⏸️ Deferred / Low Priority

### Linear Algebra

#### E-2: Use Cholesky for Symmetric Positive Definite Systems
**Priority:** Medium
**Status:** ⏸️ DEFERRED
**Reason:** Requires problem-dependent matrix type detection. Not all FEM problems produce SPD matrices.

#### E-3: Integrate MKLPardiso for Sparse Solves
**Priority:** High
**Status:** ⏸️ DEFERRED
**Reason:** Requires external dependency (Pardiso.jl or MKLSparse.jl). Consider for future integration.

#### E-4: Use Iterative Solvers for Large Systems
**Priority:** Medium
**Location:** `src/DirectXUA.jl:489`
**Status:** ⏸️ DEFERRED
**Proposed:** Use GMRES with ILU preconditioner for systems >100k DOF.

#### E-6: Optimize Eigenvalue Computation
**Priority:** Low
**Location:** `src/EigX.jl`
**Status:** ⏸️ DEFERRED
**Proposed:** Use KrylovKit's `eigsolve` for large sparse systems seeking few eigenvalues.

### FFT Optimization

#### D-4: Optimize FFT Twiddle Factor Computation
**Priority:** Medium
**Location:** `src/FFT.jl:28-38, 56-68`
**Status:** ⏸️ DEFERRED
**Issue:** Comment in line 64: "TODO use a precomputed twiddle instead".
**Proposed:** Pre-compute all twiddle factors and cache them in FFTPlan structure.

### Cache Efficiency

#### F-1: Analyze and Optimize Memory Access Patterns
**Priority:** Medium
**Status:** ⏸️ SKIPPED
**Reason:** Profiling shows no special tooling needed. Manual analysis yields limited ROI.

#### F-2: Implement Element Reordering for Cache Locality
**Priority:** Medium
**Status:** ⏸️ DEFERRED
**Reason:** Requires profiling analysis to justify implementation effort.

#### F-3: Use Array of Structures of Arrays (AoSoA)
**Priority:** Low
**Status:** ⏸️ DEFERRED
**Reason:** Major refactor. Only pursue if profiling shows significant cache issues.

#### F-4: Pre-fetch Data in Assembly
**Priority:** Low
**Status:** ⏸️ DEFERRED
**Reason:** Limited benefit on modern CPUs with hardware prefetching.

### Testing Infrastructure

#### H-2: Add Performance Regression Tests
**Priority:** High
**Status:** ✅ COMPLETED (2026-03-03)

**Implementation:**
- Created `benchmark/bench_regression.jl` - Performance regression test framework with baseline comparison
- Configurable degradation threshold (default: 20% slowdown allowed)
- Supports `--update` flag to create new baseline
- Simple text-based format (no external dependencies)
- Created GitHub Actions workflow `.github/workflows/performance.yml`:
  - Runs benchmarks on PR and compares against baseline from main
  - Posts results as PR comments with performance changes
  - Fails if performance degrades beyond threshold

**Files Created:**
- `benchmark/bench_regression.jl` - Regression test script
- `benchmark/baseline.txt` - Initial baseline performance data

#### H-3: Add Code Profiling Utilities
**Priority:** Medium
**Location:** New file `src/Profiling.jl` proposed
**Status:** ⏸️ DEFERRED
**Reason:** Not yet implemented.

#### H-4: Add Memory Allocation Tracking
**Priority:** Medium
**Status:** ⏸️ DEFERRED
**Reason:** Not yet implemented.

#### H-5: Document Performance Best Practices
**Priority:** Medium
**Status:** ⏸️ DEFERRED
**Proposed Location:** `docs/src/performance.md`

#### H-6: Add @noinline for Debug Boundaries
**Priority:** Low
**Status:** ⏸️ DEFERRED
**Reason:** Excessive inlining may bloat code and hurt instruction cache. Debug/tracing functions should not be inlined.

---

## 🔜 Future Work (Major Undertakings)

### GPU Acceleration

#### G-1: GPU-Accelerated FFT
**Priority:** High
**Location:** `src/FFT.jl`
**Status:** 🔜 FUTURE WORK
**Estimated Effort:** High (2-4 weeks)
**Proposed:** Add CUDA.jl FFT support with optional GPU path.

#### G-2: GPU-Accelerated Sparse Linear Algebra
**Priority:** High
**Location:** All solver files
**Status:** 🔜 FUTURE WORK
**Estimated Effort:** High (3-6 weeks)
**Proposed:** Use CUDA sparse solvers (CUSPARSE, CUSOLVER).

#### G-3: GPU-Accelerated Element Assembly
**Priority:** Medium
**Location:** Element assembly system
**Status:** 🔜 FUTURE WORK
**Estimated Effort:** High (4-8 weeks)
**Challenges:**
- StaticArrays work on GPU, but type parameterization needs care
- Sparse assembly requires atomic operations or specialized algorithms
- Memory transfer overhead

#### G-4: Create GPU-Compatible Element Interface
**Priority:** Low
**Location:** `src/ElementAPI.jl`
**Status:** 🔜 FUTURE WORK
**Estimated Effort:** Medium (1-2 weeks)

### Code Quality Refactoring

#### Replace Unicode Aliases with ASCII
**Priority:** High
**Location:** `src/Dialect.jl`
**Status:** 🔜 FUTURE WORK
**Estimated Effort:** Medium (1-2 weeks)

**Current:**
```julia
const ℝ = Real
const 𝕣 = Vector{Real}
const 𝕓 = Bool
const 𝕫 = Vector{Int64}
```

**Proposed:** Use standard Julia types (`Real`, `Float64`, `Int64`, `Bool`) for better portability.

**Impact:** While visually appealing, Unicode aliases make code harder to read and may cause issues in some environments.

#### Simplify @espy Macro
**Priority:** High
**Location:** `src/Espy.jl`
**Status:** 🔜 FUTURE WORK
**Estimated Effort:** High (3-4 weeks)
**Issue:** Macro is complex (400+ lines) implementing a custom DSL, creating maintenance burden.

#### Standardize Naming Conventions
**Priority:** Medium
**Status:** 🔜 FUTURE WORK
**Estimated Effort:** Medium (2-3 weeks)
**Issue:** Mixed patterns make code harder to read:
- snake_case for functions: `assemblebig!()`
- camelCase for structs: `AssemblySweepX`
- Mixed case: `∂ℝ` (starts with symbol)

#### Extract Magic Numbers to Constants
**Priority:** Low
**Status:** 🔜 FUTURE WORK
**Estimated Effort:** Low (1 week)
**Issue:** Hardcoded numerical thresholds without documentation:
```julia
const MIN_ELEMENTS_FOR_THREADING = 50
const nXnod = 2
const nUdof = 3
```

### FFT Implementation

#### Replace Custom FFT with Standard Library
**Priority:** Medium
**Status:** 🔜 FUTURE WORK
**Location:** `src/FreqXU.jl`, `src/FFT.jl`
**Estimated Effort:** Medium (2-3 weeks)
**Issue:** Custom FFT implementation (`𝔉!`, `𝔉⁻¹!`) has commented TODOs.
**Recommendation:** Consider using FFTW.jl or similar for better performance and correctness.

---

## Architecture & Code Quality

### Current Strengths

| Area | Status | Notes |
|-------|--------|-------|
| **Automatic Differentiation** | ✅ Strong | Custom AD with precedence system and StaticArrays enables efficient compilation |
| **Sparse Matrix Assembly** | ✅ Strong | Pre-computed sparsity patterns, index mapping, efficient block assembly |
| **Memory Management** | ✅ Strong | Pre-allocated buffers, in-place operations, view-based access |
| **SIMD Usage** | ✅ Fixed | Conditional outside SIMD loop for proper vectorization |
| **Multi-threading** | ✅ Recent | Thread-local assembly added (MIN_ELEMENTS_FOR_THREADING threshold) |
| **Type Stability** | 🟡 Mixed | Some issues in setup/barrier functions |
| **Multiple Solver Options** | ✅ Strong | Different solvers for different problem types (static, dynamic, optimization, frequency-domain, eigenvalue) |

### Architectural Strengths

1. **Excellent Design** - The separation of concerns (AD, assembly, solvers, elements) is well-architected

2. **Performance-Oriented** - Heavy use of StaticArrays, pre-allocated buffers, and in-place operations

3. **Extensible** - The element API allows users to add custom elements and integrate them seamlessly

### Critical Technical Debt

| Debt | Severity | Description |
|--------|------------|-------------|
| FFT implementation | Medium | Custom FFT has TODOs and may not be optimal |
| Espy macro complexity | High | Macro is hard to maintain and debug |
| Threading coverage | Low | Recent addition needs more comprehensive coverage |
| Test coverage | High | Many modules have minimal or no tests |
| Error handling | Medium | Generic error messages, no structured types |
| Documentation | Medium | Limited external documentation, no architecture diagrams |

### Potential Bugs & Reliability Issues

#### Bug: Matrix Transpose in DirectXUA
**Severity:** Medium (Correctness issue, but low impact)
**Location:** `src/DirectXUA.jl`, line 336
**Status:** ✅ INVESTIGATED (2026-03-03)

**Resolution:** The transpose implementation is mathematically correct. The forward/backward difference symmetry ensures proper index mapping between time steps.
**Validation:** Created test file `test/TestDirectXUATimeDerivativeCost.jl` that specifically tests costs on time derivatives of X and U.
**Test Result:** All 745 tests passed, including the new time derivative cost test.
**Action:** Updated TODO comment to reflect that the code has been verified.

---

## Recommendations for Next Steps

### Immediate (High Priority)

1. **Large model benchmarking** - Test with real-world models (1000+ elements with complex physics) to validate threading benefits. Create models with varied element types and physics to better understand solver performance characteristics.
2. **Profile cache performance** - Use profiling tools (Profile.jl, JET.jl) on large problems to identify actual bottlenecks, especially focusing on cache efficiency and memory access patterns.
3. **Consider iterative solvers** - For problems > 10k DOF, evaluate iterative solvers (GMRES with ILU preconditioner) as direct solvers may become expensive.
4. **Fix known issues** - Address any identified bugs or type stability issues that could affect correctness.

### Medium Term

5. **Evaluate GPU acceleration** - If CPU optimization plateaus and large-scale performance is needed for production workloads, investigate GPU acceleration for FFT and sparse linear algebra.
6. **Improve test coverage** - Add stress tests, convergence tests, and error path coverage to improve robustness.
7. **Add CI/CD pipeline** - Implement automated testing and performance regression detection for continuous integration.

### Long Term

8. **Code refactoring** - Replace Unicode aliases, simplify @espy macro, standardize naming conventions for improved maintainability.
9. **Documentation** - Add API documentation, performance guide, and architecture diagrams for better developer experience.

---

## References

### Julia Performance Resources
1. [Julia Performance Tips](https://docs.julialang.org/en/v1/manual/performance-tips/)
2. [StaticArrays.jl](https://github.com/JuliaArrays/StaticArrays.jl)
3. [LoopVectorization.jl](https://github.com/JuliaSIMD/LoopVectorization.jl)
4. [Pardiso.jl](https://github.com/JuliaSparse/Pardiso.jl)
5. [CUDA.jl](https://github.com/JuliaGPU/CUDA.jl)
6. [BenchmarkTools.jl](https://github.com/JuliaCI/BenchmarkTools.jl)
7. [Profile.jl](https://docs.julialang.org/en/v1/manual/profile/)

### FEM Optimization References
1. "High-Performance Finite Element Modeling" - T. J. R. Hughes
2. "Optimization of Sparse Matrix Assembly" - Various papers
3. "GPU Computing for Finite Elements" - GPU Gems series

---

**Document End**
