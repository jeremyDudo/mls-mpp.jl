function outer_product(a, b) # transposed, as for taichi's convention
    fin = zeros(2,2)
    fin[1] = a[1]*b[1]
    fin[2] = a[2]*b[1]
    fin[3] = a[1]*b[2]
    fin[4] = a[2]*b[2] 
    return fin
end