include("mls-mpp.jl");
using Plots


function display()
    plot = ;
end

function step!(square)
    advance()
    scatter!(square)
end

function main()

    # setup 
    square1 = add_rnd_square([0.55,0.45], 0xED553B);
    square2 = add_rnd_square([0.45,0.65], 0xF2B134);
    square3 = add_rnd_square([0.55,0.85], 0x168587);



end