include("mulMat.jl")
function polar_decomp(m) # transposed as in taichi
    x = m[1] + m[4];
    y = m[3] - m[2];
    scale = 1.0 / sqrt(x * x + y * y);
    c = x * scale;
    s = y * scale;
    R = zeros(2,2);
    R[1] = c;
    R[2] = s;
    R[3] = -s;
    R[4] = c;

    S = mulMat(m, R);

    return R, S;
end