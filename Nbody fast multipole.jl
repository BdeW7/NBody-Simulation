using Plots

mutable struct Body
    position::ComplexF64
    velocity::ComplexF64
    mass::Float64
end

mutable struct QuadNode
    center::ComplexF64
    halfwidth::Float64
    ids::Vector{Int}
    children::Vector{QuadNode}
    moments::Vector{ComplexF64}
end

const G = 1.0
const softening = 1e-3
const bucketSize = 4
const maxDepth = 20
const order = 2



function build_node(bodies, ids, center, halfwidth, depth)
    moments = zeros(ComplexF64, order + 1)
    for i in ids
        dz = bodies[i].position - center
        power = 1.0 + 0im
        for k in 0:order
            moments[k+1] += bodies[i].mass * power
            power *= dz
        end
    end
    node = QuadNode( center, halfwidth, ids, QuadNode[], moments)
    if length(ids) <= bucketSize || depth >= maxDepth
        return node
    end
    if halfwidth < 1e-12
        return node
    end
    groups = [Int[] for _ in 1:4]
    for i in ids
        z = bodies[i].position
        right = real(z) >= real(center)
        upper = imag(z) >= imag(center)
        quadrant = right ? (upper ? 4 : 2) : (upper ? 3 : 1)
        push!(groups[quadrant], i)
    end
    child_half = halfwidth / 2
    offsets = [-1-1im, 1-1im,-1+1im,1+1im]
    for q in 1:4
        if !isempty(groups[q])
            child_center = center + child_half * offsets[q]
            child = build_node(
                bodies,
                groups[q],
                child_center,
                child_half,
                depth + 1
            )
            push!(node.children, child)
        end
    end
    return node
end

function build_tree(bodies)
    ids = collect(eachindex(bodies))
    xmin = minimum(real(b.position) for b in bodies)
    xmax = maximum(real(b.position) for b in bodies)
    ymin = minimum(imag(b.position) for b in bodies)
    ymax = maximum(imag(b.position) for b in bodies)
    center = complex((xmin+xmax)/2, (ymin+ymax)/2)
    halfwidth = max((xmax-xmin)/2,(ymax-ymin)/2,1e-6) * 1.01
    return build_node(bodies, ids, center, halfwidth, 0)
end


function node_acceleration(node, bodies, target; theta=0.5)
    z = bodies[target].position
    d = z - node.center
    distance = abs(d)
    contains_target = target in node.ids
    accept = !isempty(node.children) && !contains_target && distance > 0 && node.halfwidth / distance < theta
    if accept
        a = 0.0 + 0.0im
        dc = conj(d)
        for k in 0:order
            a -= G * conj(node.moments[k+1]) /
                 dc^(k+1)
        end
        return a
    elseif isempty(node.children)
        a = 0.0 + 0.0im
        for j in node.ids
            if j == target
                continue
            end
            dz = z - bodies[j].position
            r2 = abs2(dz) + softening^2
            a -= G * bodies[j].mass * dz / r2
        end
        return a
    else
        a = 0.0 + 0.0im
        for child in node.children
            a += node_acceleration(child, bodies, target;theta=theta)
        end
        return a
    end
end


function calculate_accelerations(bodies; theta=0.5)
    tree = build_tree(bodies)
    accelerations = zeros(ComplexF64, length(bodies))
    for i in eachindex(bodies)
        accelerations[i] = node_acceleration(tree, bodies, i;theta=theta)
    end
    return accelerations
end


function step!(bodies, dt; theta=0.5)
    accelerations = calculate_accelerations(bodies; theta=theta)
    for i in eachindex(bodies)
        bodies[i].velocity += accelerations[i] * dt
        bodies[i].position += bodies[i].velocity * dt
    end
end


bodies = [Body(-0.5 + 0.0im,  0.0 + 0.7im, 10.0),Body( 0.5 + 0.0im,  0.0 - 0.7im, 10.0),Body( 0.0 + 0.5im, -0.7 + 0.0im,  1.0),Body( 0.0 - 0.5im,  0.7 + 0.0im,  1.0)]


dt = 0.001
steps = 1000


function plotter(bodies)
    x = [real(x.position) for x in bodies]
    y = [imag(y.position) for y in bodies]
    plt = scatter(
        x, y, 
        xlims=(-5,5), ylims=(-5,5), 
        markersize=1,
        legend=false
    )
end

@gif for n in 1:steps
    step!(bodies, dt; theta=0.5)
    plotter(bodies)
end every 1
