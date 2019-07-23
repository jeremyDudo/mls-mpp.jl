include("particles_to_grid.jl")
include("mod_grid_velocities.jl")
include("grid_to_particles.jl")

function advance!(party, dt)
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
    #*temp
    F = [1 0; 0 1];  # Deformation tensor
    #*temp
    C = [0 0; 0 0];  # Cauchy tensor
    Jp = 1;          # Jacobian determinant (scalar)
    # Reset grid
    grid = zeros((n+1)^2,3)

    # 1. Particles to grid
    particles_to_grid!(party, n, grid, dx, inv_dx, dt, vol, hardening, particle_mass, μ₀, λ₀)

    # Modify grid velocities to respect boundaries
    boundary = 0.05;
    mod_grid_velocities!(n, grid, dt, boundary)

    # 2. Grid to particle
    grid_to_particles!(party, n, grid, inv_dx, plastic)

end