using Profile
using Muscade
using Muscade.Toolbox
using StaticArrays

@functor with() load(t) = 1e5 * t

function run_benchmark()
    model = Model(:BenchModel)

    n_ele = 1000
    L = 10.0
    nnodes = n_ele + 1

    nodeCoord = hcat(collect(range(0, L, length=nnodes)), zeros(nnodes), zeros(nnodes))
    nodes = addnode!(model, nodeCoord)

    mat = BeamCrossSection(EA=1e9, EI₂=1e8, EI₃=1e8, GJ=1e7, μ=100.0, ι₁=10.0, Ca₂=0.0, Cq₂=0.0, Ca₃=0.0, Cq₃=0.0)

    mesh = hcat(nodes[1:n_ele], nodes[2:n_ele+1])
    addelement!(model, EulerBeam3D, mesh; mat=mat, orient2=SVector(0., 1., 0.))

    addelement!(model, Hold, [nodes[1]]; field=:t1)
    addelement!(model, Hold, [nodes[1]]; field=:t2)
    addelement!(model, Hold, [nodes[1]]; field=:t3)
    addelement!(model, Hold, [nodes[1]]; field=:r1)
    addelement!(model, Hold, [nodes[1]]; field=:r2)
    addelement!(model, Hold, [nodes[1]]; field=:r3)

    # Add a load to bend it
    addelement!(model, DofLoad, [nodes[end]]; field=:t2, value=load)

    initialstate = initialize!(model)
    time_steps = 0.0:0.1:1.0

    # Run SweepX{0} - Newton Raphson static
    solve(SweepX{0}; initialstate, time=time_steps, verbose=false, maxΔx=1e-5, maxiter=50)
end

println("Compiling benchmark...")
run_benchmark()

println("Profiling benchmark...")
Profile.clear()
@profile run_benchmark()

open("profile_results.txt", "w") do io
    Profile.print(io, format=:flat, sortedby=:count, mincount=1)
end

println("Done profiling.")
