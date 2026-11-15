using Plots

# Get the normalised direction vector
function norm(vec)
    S = sqrt(vec[1]^2 + vec[2]^2 + vec[3]^2) #Magnitude of vector
    return vec / S
end

# Find Force on object one fromm object two given positions and mass (G is gravity variable).
function FindOneForce(pos1,pos2, mas1,mas2; G=1)

    r = (pos1[1]-pos2[1])^2 + (pos1[2]-pos2[2])^2 + (pos1[3]-pos2[3])^2
    #Square of distance: the plain distance is not needed, so don't bother with the root

    F = G*mas1*mas2/(r)
    #Magnitude of Force on obj 1

    dir = norm(pos2 .- pos1)
    #Normalised direction vector for the force

    return dir * F
end

#Find all forces between a set of objects, given position and masses.
function FindForces(pos,mas)
    n=length(pos)
    #Number of objects

    forces = [[0.0,0.0,0.0] for i in 1:n]
    #Initialised list of forces, starting at zero so we can add up the force vectors on the objects

    #Since the force on i from j is the negative of the force on j from i, we only need to calculate for i<j.
    for i in 1:(n-1) 
        for j in (i+1):(n)
            F = FindOneForce(pos[i],pos[j],mas[i],mas[j])
            forces[i] = forces[i] + F
            forces[j] = forces[j] - F
        end
    end

    return forces
end


#Datastructure for storing the information about the objects
mutable struct Nbody
    bodies::Int64
    positions::Vector{Vector{Float64}}
    masses::Vector{Float64}
    velocities::Vector{Vector{Float64}}

    forces::Vector{Vector{Float64}}
end

#Initialise a random Nbody system
function InitNbody(n; positions=[collect(t) for t in zip(rand(n),rand(n),rand(n))], masses=rand(n), velocities=[collect(t) for t in zip(rand(n),rand(n),rand(n))])
    forces = FindForces(positions, masses)
    return Nbody(n,positions,masses,velocities, forces)
end

#Move forward by a tick in the simulation
function step!(body::Nbody)
    n= body.bodies

    #update velocities by adding the acceleration (F=ma => a=F/m)
    #then update positions by these velocities
    for i in 1:n
        body.velocities[i] .+= (body.forces[i] / body.masses[i])
        body.positions[i] .+= body.velocities[i]
    end

    body.forces=FindForces(body.positions, body.masses)
    #update the forces of the system
end

#plot the system
function plotter(body::Nbody)
    v = body.positions
    x = getindex.(v,1)
    y = getindex.(v,2)
    z = getindex.(v,3)
    plt = scatter(
        x, y, z,
        xlims=(-1,3), ylims=(-1,3), zlims=(-1,3),
        markersize=5,
        legend=false
    )
end



x=0.0227
n=5
p = [[2.0,1.5,1.0],[1.0,1.5,1.0], [1.0,1.5,2.0], [2.0,1.5,2.0],[1.5,1.5,1.5]]
m = [x,x,x,x,1.0]
v = [[0.1,0.0,0.1],[0.1,0.0,-0.1],[-0.1,0.0,-0.1],[-0.1,0.0,0.1],[0.0,0.0,0.0]]


example = InitNbody(n; positions=p, masses=m, velocities=v)

@gif for t=1:20
    step!(example)
    plotter(example)
end every 1
