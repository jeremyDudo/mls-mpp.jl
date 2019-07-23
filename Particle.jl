mutable struct Particle
    x # position
    v # velocity
    F # Deformation tensor
    C # Cauchy tensor
    Jp # Jacobian determinant (scalar)
    c # color (int)
end