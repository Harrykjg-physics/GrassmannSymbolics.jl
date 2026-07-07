import Pkg

Pkg.activate(@__DIR__)
include(joinpath(
    @__DIR__,
    string(:examples),
    string(:Single_Flavor_Gross_Neveu_Wilson, Char(46), :jl),
))

model = derive_single_flavor_gnw_tensor()
D, g2, z, s = model.parameters
reference = published_gnw_tensor(model.legs, D, g2, z, s)
println((length(model.tensor.terms), length(reference.terms)))

for mask in 0:255
    bits = gnw_occupation_bits(mask)
    occupied = gnw_occupied_legs(model.legs, bits)
    actual = get_coeff(model.tensor, occupied)
    expected = published_gnw_coefficient(bits, D, g2, z, s)
    if !symbolically_equal(GrassmannExpr(actual), GrassmannExpr(expected))
        difference = simplify_coefficients(GrassmannExpr(actual - expected))
        println((mask, bits, actual, expected, difference))
    end
end
