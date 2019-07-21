include("mulMat.jl")
function polar_decomp(m) # transposed as in taichi
    x = m[0] + m[3];
    y = m[2] - m[1];
    scale = 1.0 / sqrt(x * x + y * y);
    c = x * scale;
    s = y * scale;
    R = [];
    R[0] = c;
    R[1] = s;
    R[2] = -s;
    R[3] = c;
    S = mulMat(m, R);

    return R, S;
end