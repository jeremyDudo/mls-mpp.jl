
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