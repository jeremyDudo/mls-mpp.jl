include("outer_product.jl")
include("gridIndex.jl")
include("mulMat.jl")
include("polar_decomp.jl")
include("mulMatVec.jl")

using LinearAlgebra, Random;

function particles_to_grid!(party, n, grid, dx, inv_dx, dt, vol, hardening, particle_mass, μ₀, λ₀)
    for p in party 
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
        
        stress =  k1 .* ( mulMat((p.F'.- r), p.F) .* 2 .* μ .+ [k2 0; 0 k2] );
        
        affine =  stress .+ (particle_mass .* p.C);
        
        mv = [p.v[1]*particle_mass, p.v[2]*particle_mass, particle_mass]; # translational momentum
        for i = 1:3, j = 1:3
            # scatter to grid
            dpos =    [(i-fx[1])*dx, (j-fx[2])*dx];
            
            # ii =      gridIndex(base_coord[1] + i, base_coord[2] + j, n);
            
            weight =  w[i][1] * w[j][2];
            grid[i,j] = (grid[i,j] .+ ((mv .+ [mulMatVec(affine, dpos)[1], mulMatVec(affine, dpos)[2], 0]) .* weight))[1];
        end
    end
end