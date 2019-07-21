function mulMatVec(a, b) # transposed, as for taichi's convention
    return [
        a[1]*b[1]+a[3]*b[2],
        a[2]*b[1]+a[4]*b[2]
    ];
end