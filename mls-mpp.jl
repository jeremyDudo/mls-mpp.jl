using LinearAlgebra, Random;
include("outer_product.jl")
include("mulMat.jl")
include("polar_decomp.jl")
include("mulMatVec.jl")
# grid resolution
n = 80;

widow_size = 800;

dt = 1e-4;
frame_dt = 1e-3;
dx = 1.0 / n;
inv_dx = 1.0 / dx;

particle_mass = 1.0;
vol = 1.0;

hardening = 10.0;
E = 1e4;
ν = 0.2;
μ₀ = E / (2 * (1 + ν));
λ₀ = E * ν / ((1+ν) * (1 - 2 * ν));
plastic = 1;

v = [0, 0];       # velocity
F = [1, 0, 0, 1];  # Deformation tensor
C = [0, 0, 0, 0];  # Cauchy tensor
Jp = 1;          # Jacobian determinant (scalar)

mutable struct Particle
    x # position
    v # velocity
    F # Deformation tensor
    C # Cauchy tensor
    Jp # Jacobian determinant (scalar)
    c # color (int)
end

particles = [];
grid = [];    # velocity + mass, node_res = cell_res + 1

gridIndex(i, j) = i + (n+1)*j;

function advance(dt) 
    # Reset grid
    grid = zeros((n+1)^2,3)

    # 1. Particles to grid
    for p in particles 
        base_coord = floor.( Int, (inv_dx .* p.x .- 0.5) ); # element-wise floor
        fx = inv_dx .* p.x .- base_coord; # base position in grid units

        # Quadratic kernels  [http:#mpm.graphics   Eqn. 123, with x=fx, fx-1,fx-2]
        w = [
            (0.5 .* map(x -> x^2, (1.5 .- fx))),
            (0.75 .- map(x -> x^2, (fx .- 1.0))),
            (0.5 .* map(x -> x^2, (fx .- 0.5)))
        ];

        # Snow-like hardening
        e = exp(hardening * (1.0 - p.Jp));
        μ = μ₀ * e;
        λ = λ₀ * e;

        # Cauchy stress times dt and inv_dx
        # original taichi: stress = -4*inv_dx*inv_dx*dt*vol*( 2*mu*(p.F-r)*transposed(p.F) + lambda*(J-1)*J )
        # (in taichi matrices are coded transposed)
        J =       det(p.F);           # Current volume
        r, s =    polar_decomp(p.F);  # Polar decomp. for fixed corotated model
        k1 =      -4*inv_dx*inv_dx*dt*vol;
        k2 =      λ*(J-1)*J;
        stress =  map(x -> x*k1, ( map(x-> x*2*μ, mulMat((p.F'.- r), p.F)) .+ [k2, 0, 0, k2] ) );
        affine =  stress .+ map(x -> x*particle_mass, p.C);

        mv = [p.v[1]*particle_mass, p.v[2]*particle_mass, particle_mass]; # translational momentum
        for i = 1:3, j = 1:3
            # scatter to grid
            dpos =    [(i-fx[1])*dx, (j-fx[2])*dx];
            ii =      gridIndex(base_coord[1] + i, base_coord[2] + j);
            weight =  w[i][1] * w[j][2];
            grid[ii] =      (grid[ii] .+ ((mv .+ [mulMatVec(affine, dpos)[1], mulMatVec(affine, dpos)[2], 0]) .* weight));
        end
    end

    # Modify grid velocities to respect boundaries
    boundary = 0.05;
    for i = 1:n, j = 1:n
        # for all grid nodes
          ii = gridIndex(i, j);
        if (grid[ii][2] > 0)  # no need for epsilon here
            grid[ii] =  map(x -> x/grid[ii][3], grid[ii]); # normalize by mass
            grid[ii] =  grid[ii] .+ [0, -200*dt, 0]; # add gravity
              x =   i/n;
              y =   j/n; # boundary thickness, node coord

            # stick
            if (x < boundary|| x > 1-boundary || y > 1-boundary) 
                grid[ii] = [0, 0, 0];
            end

            # separate
            if (y < boundary) 
                grid[ii][2] = max(0.0, grid[ii][2]);
            end
        end
    end

    # 2. Grid to particle
    for p in particles
        base_coord= floor( Int, (inv_dx .* p.x .- 0.5));# element-wise floor
        fx = ((p.x .* inv_dx) .- base_coord); # base position in grid units
        w = [
            (0.5 .* map(x -> x^2, (1.5 .- fx))),
            (0.75 .- map(x -> x^2, (fx .- 1.0))),
            (0.5 .* map(x -> x^2, (fx .- 0.5)))
        ];
        p.C = [0,0, 0,0];
        p.v = [0, 0];
        for i = 1:3, j = 1:3 
            dpos =    [i, j] .- fx;
            ii =      gridIndex(base_coord[1] + i, base_coord[2] + j);
            weight =  w[i][1] * w[j][2];
            p.v =     (p.v .+ (grid[ii] .* weight)); # velocity
            p.C =     map(x -> x*4*inv_dx, p.C .+ outer_product((grid[ii] .* weight), dpos)); # APIC (affine particle-in-cell); p.C is the affine momentum
        end

        # advection
        p.x = (p.x .+ (p.v .* dt));

        # MLS-MPM F-update
        # original taichi: F = (Mat(1) + dt * p.C) * p.F
        F = mulMat(p.F, ([1,0, 0,1] .+ dt .* p.C));

        # Snow-like plasticity
        M = svd(F);
        svd_u, sig, svd_v = M.U, M.S, M.V;
        for i = 1:2*plastic  
            sig[i+2*i] = clamp(sig[i+2*i], 1.0 - 2.5e-2, 1.0 + 7.5e-3);
        end
          oldJ = det(F);
        # original taichi: F = svd_u * sig * transposed(svd_v)
        F = mulMat(mulMat(svd_u, sig), svd_v');
          Jp_new = clamp(p.Jp * oldJ / det(F), 0.6, 20.0);
        p.Jp = Jp_new;
        p.F = F;
    end
end

function add_rnd_square(center, c) 
    for i = 1:1000 
        # Randomly sample 1000 particles in the square
        push!(particles, (Particle(([(rand()*2-1)*0.08, (rand()*2-1)*0.08] .+ center), v, F, C, Jp, c)));
    end
end

function add_rnd_square2(center, c)
    party = []
    for i = 1:1000 
        # Randomly sample 1000 particles in the square
        push!(party, (Particle(([(rand()*2-1)*0.08, (rand()*2-1)*0.08] .+ center), v, F, C, Jp, c)));
    end
    party
end