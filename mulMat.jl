function mulMat(a, b) 
    return [
        a[1]*b[1]+a[2]*b[3],
        a[1]*b[2]+a[2]*b[4],
        a[3]*b[1]+a[4]*b[3],
        a[3]*b[2]+a[4]*b[4]
    ];
end