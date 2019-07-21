using LinearAlgebra;
# grid resolution
const n = 80;

const widow_size = 800;

const dt = 1e-4;
const frame_dt = 1e-3;
const dx = 1.0 / n;
const inv_dx = 1.0 / dx;

const particle_mass = 1.0;
const vol = 1.0;

const hardening = 10.0;
const E = 1e4;
const ν = 0.2;
const μ₀ = E / (2 * (1 + ν));
const λ₀ = E * ν / ((1+ν) * (1 - 2 * ν));

v = [0, 0];       # velocity
F = [1, 0, 0, 1];  # Deformation tensor
C = [0, 0, 0, 0];  # Cauchy tensor
Jp = 1;          # Jacobian determinant (scalar)

mutable struct Particle
    x, # position
    v, # velocity
    F, # Deformation tensor
    C, # Cauchy tensor
    J, # Jacobian determinant (scalar)
    c # color (int)
end

const particles = [];
const grid = [];    # velocity + mass, node_res = cell_res + 1

gridIndex(i, j) = i + (n+1)*j;

function advance(dt) 
    # Reset grid
    grid = zeros((n+1)^2,3)


    # 1. Particles to grid
    for p in particles 
        const base_coord = map(x -> Int(round(x)), (inv_dx .* p.x .- 0.5) ); # element-wise floor
        const fx = inv_dx .* p.x .- base_coord; # base position in grid units

        # Quadratic kernels  [http:#mpm.graphics   Eqn. 123, with x=fx, fx-1,fx-2]
        const w = [
            map(x -> x^2, (0.5 .* (1.5 .- fx))),
            map(x -> x^2, (0.75 .- (fx .- 1.0))),
            map(x -> x^2, (0.5 .* (fx .- 0.5)))
        ];

        # Snow-like hardening
        const e = exp(hardening * (1.0 - p.Jp));
        const μ = μ₀ * e;
        const λ = λ₀ * e;

        # Cauchy stress times dt and inv_dx
        # original taichi: stress = -4*inv_dx*inv_dx*dt*vol*( 2*mu*(p.F-r)*transposed(p.F) + lambda*(J-1)*J )
        # (in taichi matrices are coded transposed)
        const J = det(p.F); # Current volume
        const {R:r, S:s} = polar_decomp(p.F); # Polar decomp. for fixed corotated model
        const k1 = -4*inv_dx*inv_dx*dt*vol;
        const k2 = λ*(J-1)*J;
        const stress = map(x -> x*k1, ( map(x-> x*2*μ, mulMat((p.F'.- r),p.F)) .+ [k2,0,0,k2] ) );
        const affine = stress .+ map(x -> x*particle_mass, p.C);

        const mv = [p.v[1]*particle_mass, p.v[2]*particle_mass, particle_mass]; # translational momentum
        for i = 1:3, j = 1:3
            # scatter to grid
            const dpos = [(i-fx[1])*dx, (j-fx[2])*dx];
            const ii = gridIndex(base_coord[1] + i, base_coord[2] + j);
            const weight = w[i][1] * w[j][2];
            grid[ii] = (grid[ii] .+ ((mv .+ [mulMatVec(affine, dpos)[1], mulMatVec(affine, dpos)[2],0]) .* weight));
        end
    end

    # Modify grid velocities to respect boundaries
    const boundary = 0.05;
    for i = 1:n, j = 1:n
        # for all grid nodes
        const ii = gridIndex(i, j);
        if (grid[ii][2] > 0)  # no need for epsilon here
            grid[ii] = map(x -> x/grid[ii][3], grid[ii]); # normalize by mass
            grid[ii] = grid[ii] .+ [0, -200*dt, 0]; # add gravity
            const x = i/n;
            const y = j/n; # boundary thickness, node coord

            # stick
            if (x < boundary||x > 1-boundary||y > 1-boundary) 
                grid[ii]=[0,0,0];
            end

            # separate
            if (y < boundary) 
                grid[ii][2] = max(0.0, grid[ii][2]);
            end
        end
    end

    # 2. Grid to particle
    for p in particles
        const base_coord=map(x -> Int(round(x)), (map(x -> x*inv_dx, p.x) .- 0.5));# element-wise floor
        const fx = ((p.x .* inv_dx) .- base_coord); # base position in grid units
        const w = [
            had2D([0.5, 0.5], sub2D([1.5, 1.5], fx).map(o=>o*o)),
            sub2D([0.75, 0.75], sub2D(fx, [1.0, 1.0]).map(o=>o*o)),
            had2D([0.5, 0.5], sub2D(fx, [0.5,0.5]).map(o=>o*o))
        ];
        p.C = [0,0, 0,0];
        p.v = [0, 0];
        for i = 1:3, j = 1:3 
            const dpos = sub2D([i, j], fx);
            const ii = gridIndex(base_coord[0] + i, base_coord[1] + j);
            const weight = w[i][0] * w[j][1];
            p.v = add2D(p.v, sca2D(grid[ii], weight)); # velocity
            p.C = addMat(p.C, outer_product(sca2D(grid[ii],weight), dpos).map(o=>o*4*inv_dx)); # APIC (affine particle-in-cell); p.C is the affine momentum
        end

        # advection
        p.x = add2D(p.x, sca2D(p.v, dt));

        # MLS-MPM F-update
        # original taichi: F = (Mat(1) + dt * p.C) * p.F
        F = mulMat(p.F, addMat([1,0, 0,1], p.C.map(o=>o*dt)));

        # Snow-like plasticity
        {U:svd_u, sig:sig, V:svd_v} = svd(F);
        for i = 1:2*plastic  
            sig[i+2*i] = clamp(sig[i+2*i], 1.0 - 2.5e-2, 1.0 + 7.5e-3);
        end
        const oldJ = determinant(F);
        # original taichi: F = svd_u * sig * transposed(svd_v)
        F = mulMat(mulMat(svd_u, sig), transposed(svd_v));
        const Jp_new = clamp(p.Jp * oldJ / determinant(F), 0.6, 20.0);
        p.Jp = Jp_new;
        p.F = F;
        end
    end