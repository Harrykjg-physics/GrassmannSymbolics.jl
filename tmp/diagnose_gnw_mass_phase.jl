import Pkg

Pkg.activate(@__DIR__)
include(joinpath(
    @__DIR__,
    string(:examples),
    string(:Single_Flavor_Gross_Neveu_Wilson, Char(46), :jl),
))

model = derive_single_flavor_gnw_tensor()
D, g2, z, s = model.parameters
for mask in 0:255
    bits = gnw_occupation_bits(mask)
    i1, j1, i2, j2, j1p, i1p, j2p, i2p = bits
    bar_degree = i1 + i2 + j1p + j2p
    unbar_degree = j1 + j2 + i1p + i2p
    if bar_degree == 1 && unbar_degree == 1
        ordering_phase = (
            i1 * (i2p + i1p + j2 + j1) +
            i2 * (i2p + i1p + j2) +
            j1p * (i2p + i1p) +
            j2p * i2p + i1p + i2p
        )
        chemical_power = i2 + i2p - j2 - j2p
        prefactor = (-1)^ordering_phase * s^2 * z^chemical_power
        actual = get_coeff(model.tensor, gnw_occupied_legs(model.legs, bits))
        normalized = Symbolics.simplify(actual / (prefactor * D))
        println((bits, normalized))
    end
end
