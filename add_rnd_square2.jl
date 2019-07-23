using Random;
include("Particle.jl")

function add_rnd_square2(center, c)
    party = []
    for i = 1:1000 
        # Randomly sample 1000 particles in the square
        push!(party, (Particle(([(rand()*2-1)*0.08, (rand()*2-1)*0.08] .+ center), v, F, C, Jp, c)));
    end
    party
end