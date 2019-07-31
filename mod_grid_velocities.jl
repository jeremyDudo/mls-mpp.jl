include("gridIndex.jl")

function mod_grid_velocities!(n, grid, dt, boundary)
    for i = 1:n, j = 1:n
        # for all grid nodes
        ii = gridIndex(i, j, n);
        if (grid[ii][3] > 0)  # no need for epsilon here
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
end