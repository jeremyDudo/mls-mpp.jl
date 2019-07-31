include("outer_product.jl")
include("gridIndex.jl")
include("mulMat.jl")

using LinearAlgebra;

function grid_to_particles!(party, n, grid, inv_dx, plastic, dt)
    for p in party
        base_coord= floor.( Int, (inv_dx .* p.x .- 0.5));# element-wise floor
        fx = ((p.x .* inv_dx) .- base_coord); # base position in grid units
        w = [
            (0.5 .* map(x -> x^2, (1.5 .- fx))),
            (0.75 .- map(x -> x^2, (fx .- 1.0))),
            (0.5 .* map(x -> x^2, (fx .- 0.5)))
        ];
        p.C = [0 0; 0 0];
        p.v = [0, 0];
        for i = 1:3, j = 1:3 
            dpos =    [i, j] .- fx;
            ii =      gridIndex(base_coord[1] + i, base_coord[2] + j, n);
            weight =  w[i][1] * w[j][2];
            println(p.C)
            println(grid[ii])
            println(weight)
            println(outer_product((grid[ii] .* weight), dpos))
            p.v =     (p.v .+ (grid[ii] .* weight)); # velocity

            p.C =     map(x -> x*4*inv_dx, p.C .+ outer_product((grid[ii] .* weight), dpos)); # APIC (affine particle-in-cell); p.C is the affine momentum
        end

        # advection
        p.x = (p.x .+ (p.v .* dt));

        # MLS-MPM F-update
        # original taichi: F = (Mat(1) + dt * p.C) * p.F
        F = mulMat(p.F, ([1 0; 0 1] .+ dt .* p.C));

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