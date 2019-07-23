include("mls-mpp.jl");
using Plots


function main()

    """
    I fear that advancing each square separately might be a bad idea, but that can be an easy fix if I can animate anything
    """
    # setup 

    dt = 0.1

    square1 = add_rnd_square2([0.55,0.45], 0xED553B);
    square2 = add_rnd_square2([0.45,0.65], 0xF2B134);
    square3 = add_rnd_square2([0.55,0.85], 0x168587);

    @gif for i = 1:dt:1500
        advance!(square1, dt)
        advance!(square2, dt)
        advance!(square3, dt)

        scatter!((square1, square2, square3)...)
    end every 10
end

main()