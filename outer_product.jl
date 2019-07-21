function outer_product(a, b) # transposed, as for taichi's convention
    return [
        a[1]*b[1],a[2]*b[1],
        a[1]*b[2],a[2]*b[2]
    ];
end