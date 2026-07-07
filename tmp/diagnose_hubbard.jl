using Pkg

slash = Char(47)
dot = Char(46)
root = @__DIR__() * string(slash, dot, dot)
Pkg.activate(root)

include(joinpath(root, string(:examples), string(1, :D_Hubbard, dot, :jl)))

model = derive_hubbard_tensor()
p = model.parameters

for mask in 0:4095
    leg_bits = hubbard_leg_bits(mask)
    bits = hubbard_bits_from_leg_bits(leg_bits)
    occupied = GrassmannGenerator[
        leg for (leg, occupation) in zip(model.legs, leg_bits) if occupation == 1
    ]
    actual = get_coeff(model.tensor, occupied)
    expected = published_hubbard_coefficient(bits, p.q, p.h, p.Ue, model.named_legs)
    if !symbolically_equal(GrassmannExpr(actual), GrassmannExpr(expected))
        println(string(:mismatch))
        println(string(:mask, :_, mask))
        println(leg_bits)
        println(bits)
        println(string(:args, :_, hubbard_delta_arguments(bits)))
        println(string(:spatial_degree, :_, hubbard_spatial_degree(bits)))
        println(string(:forward_degree, :_, hubbard_forward_spatial_degree(bits)))
        println(string(:Rsign, :_, hubbard_R_sign(bits, model.named_legs)))
        println(string(:actual, :_, actual))
        println(string(:expected, :_, expected))
        println(string(:sum, :_, Symbolics.simplify(actual + expected; expand=true)))
        println(string(:difference, :_, Symbolics.simplify(actual - expected; expand=true)))
        exit()
    end
end

println(string(:no_mismatch))
