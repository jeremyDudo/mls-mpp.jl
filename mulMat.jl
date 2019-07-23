function mulMat(a, b) 
    # println(a)
    # println(b)
    fin = zeros(2,2)
    fin[1] = a[1]*b[1]+a[2]*b[3]
    fin[2] = a[1]*b[2]+a[2]*b[4]
    fin[3] = a[3]*b[1]+a[4]*b[3]
    fin[4] = a[3]*b[2]+a[4]*b[4]
    return fin
end