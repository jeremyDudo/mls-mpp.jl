function advance(dt) 
    # Reset grid
    for(let i = 0; i < (n+1)*(n+1); i++) 
        grid[i] = [0,0,0];
    end

    # 1. Particles to grid
    for (let p of particles) 
        const base_coord=sub2D(sca2D(p.x, inv_dx), [0.5,0.5]).map((o)=>parseInt(o)); # element-wise floor
        const fx = sub2D(sca2D(p.x, inv_dx), base_coord); # base position in grid units

        # Quadratic kernels  [http:#mpm.graphics   Eqn. 123, with x=fx, fx-1,fx-2]
        const w = [
            had2D([0.5, 0.5], sub2D([1.5, 1.5], fx).map(o=>o*o)),
            sub2D([0.75, 0.75], sub2D(fx, [1.0, 1.0]).map(o=>o*o)),
            had2D([0.5, 0.5], sub2D(fx, [0.5, 0.5]).map(o=>o*o))
        ];

        # Snow-like hardening
        const e = Math.exp(hardening * (1.0 - p.Jp));
        const mu=mu_0*e;
        const lambda=lambda_0*e;

        # Cauchy stress times dt and inv_dx
        # original taichi: stress = -4*inv_dx*inv_dx*dt*vol*( 2*mu*(p.F-r)*transposed(p.F) + lambda*(J-1)*J )
        # (in taichi matrices are coded transposed)
        const J = determinant(p.F); # Current volume
        const {R:r, S:s} = polar_decomp(p.F); # Polar decomp. for fixed corotated model
        const k1 = -4*inv_dx*inv_dx*dt*vol;
        const k2 = lambda*(J-1)*J;
        const stress = addMat( mulMat(subMat(transposed(p.F),r),p.F).map(o=>o*2*mu), [k2,0,0,k2] ).map(o=>o*k1);
        const affine = addMat(stress, p.C.map(o=>o*particle_mass));

        const mv = [p.v[0]*particle_mass, p.v[1]*particle_mass, particle_mass]; # translational momentum
        for (let i = 0; i < 3; i++) 
            for (let j = 0; j < 3; j++)  # scatter to grid
                const dpos = [(i-fx[0])*dx, (j-fx[1])*dx];
                const ii = gridIndex(base_coord[0] + i, base_coord[1] + j);
                const weight = w[i][0] * w[j][1];
                grid[ii] = add3D(grid[ii], sca3D(add3D(mv, [...mulMatVec(affine, dpos),0]), weight));
            end
        end
    end

    # Modify grid velocities to respect boundaries
    const boundary = 0.05;
    for(let i = 0; i <= n; i++) 
        for(let j = 0; j <= n; j++)  # for all grid nodes
            const ii = gridIndex(i, j);
            if (grid[ii][2] > 0)  # no need for epsilon here
                grid[ii] = grid[ii].map(o=>o/grid[ii][2]); # normalize by mass
                grid[ii] = add3D(grid[ii], [0,-200*dt,0]); # add gravity
                const x = i/n;
                const y = j/n; # boundary thickness, node coord

                # stick
                if (x < boundary||x > 1-boundary||y > 1-boundary) 
                    grid[ii]=[0,0,0];
                end

                # separate
                if (y < boundary) 
                    grid[ii][1] = Math.max(0.0, grid[ii][1]);
                end
            end
        end
    end

    # 2. Grid to particle
    for (let p of particles) 
        const base_coord=sub2D(p.x.map(o=>o*inv_dx),[0.5,0.5]).map(o=>parseInt(o));# element-wise floor
        const fx = sub2D(sca2D(p.x, inv_dx), base_coord); # base position in grid units
        const w = [
            had2D([0.5, 0.5], sub2D([1.5, 1.5], fx).map(o=>o*o)),
            sub2D([0.75, 0.75], sub2D(fx, [1.0, 1.0]).map(o=>o*o)),
            had2D([0.5, 0.5], sub2D(fx, [0.5,0.5]).map(o=>o*o))
        ];
        p.C = [0,0, 0,0];
        p.v = [0, 0];
        for (let i = 0; i < 3; i++) 
            for (let j = 0; j < 3; j++) 
                const dpos = sub2D([i, j], fx);
                const ii = gridIndex(base_coord[0] + i, base_coord[1] + j);
                const weight = w[i][0] * w[j][1];
                p.v = add2D(p.v, sca2D(grid[ii], weight)); # velocity
                p.C = addMat(p.C, outer_product(sca2D(grid[ii],weight), dpos).map(o=>o*4*inv_dx)); # APIC (affine particle-in-cell); p.C is the affine momentum
            end
        end

        # advection
        p.x = add2D(p.x, sca2D(p.v, dt));

        # MLS-MPM F-update
        # original taichi: F = (Mat(1) + dt * p.C) * p.F
        let F = mulMat(p.F, addMat([1,0, 0,1], p.C.map(o=>o*dt)));

        # Snow-like plasticity
        let {U:svd_u, sig:sig, V:svd_v} = svd(F);
        for (let i = 0; i < 2 * plastic; i++) 
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