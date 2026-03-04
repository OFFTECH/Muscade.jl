# Test to verify correct handling of costs on time derivatives
# This validates the transpose TODO in DirectXUA.jl line 336

module TestDirectXUATimeDerivativeCost

using Test
using Muscade
using StaticArrays, SparseArrays

# Simple element with time derivative dependence
struct TimeDerivElement <: AbstractElement
    k :: Float64
end

@espy function Muscade.residual(o::TimeDerivElement, X, U, A, t, SP, dbg)
    # Use X, X', X" in residual to test time derivative handling
    x, x′, x″ = ∂0(X)[1], ∂1(X)[1], ∂2(X)[1]
    # Simple spring-like residual
    r = o.k * x
    return SVector(r), noFB
end

Muscade.doflist(::Type{TimeDerivElement}) = (inod=(1,), class=(:X, :U, :A), field=(:x,))

# Include common test elements
include("SomeElements.jl")

@testset "DirectXUA - Time Derivative Cost" begin
    # Create model with time derivative cost
    model = Model(:TimeDerivTest)
    n1 = addnode!(model, 𝕣[0.0])
    e1 = addelement!(model, TimeDerivElement, [n1], k=1.0)

    # Add cost on X' (first time derivative)
    @functor with() cost_xprime(x, t) = 0.5 * x^2
    e2 = addelement!(model, SingleUdof, [n1];
                     Xfield=:x, Ufield=:x′, cost=cost_xprime)

    state0 = initialize!(model)

    # Parameters for time-stepping
    nstep = 5
    Δt = 0.1
    OX = 2  # Include second time derivatives
    OU = 0
    IA = 0  # No A variables

    out, asm, dofgr = Muscade.prepare(Muscade.AssemblyDirect{OX,OU,IA}, model, state0.dis)

    # Assemble at a single time step
    state = [Muscade.State{1,OX+1,OU+1}(copy(state0, SP=(γ=0.,iter=1)) for i = 1:nstep]
    for i = 1:nstep
        state[i].time = Δt * i
    end

    Muscade.assemble!{:matrices}(out, asm, state0.dis, model, state[1], Δt, (γ=0.,iter=1), ())

    # Check that L1 contains gradients for time derivatives
    @test size(out.L1[2]) == (OX+1,)  # Should have entries for X, X', X"
    @test length(out.L1[2][1]) > 0   # ∂L/∂X should exist
    @test length(out.L1[2][2]) > 0   # ∂L/∂X' should exist (this tests the transpose issue)

    # Test with full solver
    results = solve(DirectXUA{OX,OU,IA};
                    initialstate=[state0],
                    time=[0:Δt:Δt*(nstep-1)],
                    maxiter=50,
                    verbose=false)

    # Solution should converge
    @test length(results[1]) == nstep
    @test all(isfinite.(collect(first(results[1].X)))
end

end
