using GrassmannSymbolics
using Symbolics
using Test

# Direct local integration of equations (13)-(23) in the reference note.
# D=m+2r, g2=g^2, z=exp(mu/2), and s=1/sqrt(2).
function derive_single_flavor_gnw_tensor()
    @variables D g2 z s

    ψ1, ψ2 = grassmann(:ψ, 1:2)
    ψbar1, ψbar2 = dual.((ψ1, ψ2))
    η1, ξ1, η2, ξ2 = (
        grassmann(:η1), grassmann(:ξ1), grassmann(:η2), grassmann(:ξ2)
    )
    ηbar1, ξbar1, ηbar2, ξbar2 = dual.((η1, ξ1, η2, ξ2))

    density1 = ψbar1 * ψ1
    density2 = ψbar2 * ψ2
    onsite_action = D * (density1 + density2) - 2g2 * density1 * density2

    bond_exponents = [
        s * (ψbar1 - ψbar2) * η1,
        s * ηbar1 * (ψ1 - ψ2),
        s * z * (ψbar1 - im * ψbar2) * η2,
        s * z * ηbar2 * (ψ1 + im * ψ2),
        s * (ψbar1 + ψbar2) * ξbar1,
        -s * ξ1 * (ψ1 + ψ2),
        (s / z) * (ψbar1 + im * ψbar2) * ξbar2,
        -(s / z) * ξ2 * (ψ1 - im * ψ2),
    ]
    tensor = local_grassmann_tensor(
        onsite_action,
        bond_exponents,
        [ψ1, ψbar1, ψ2, ψbar2],
    )

    # Flat form of (η1 ξ1)(η2 ξ2)(ξbar1 ηbar1)(ξbar2 ηbar2).
    legs = (η1, ξ1, η2, ξ2, ξbar1, ηbar1, ξbar2, ηbar2)
    parameters = (; D, g2, z, s)
    return (; tensor, legs, parameters)
end

function gnw_Abar(bits)
    bits == (1, 1, 0, 0) && return -1 + im
    bits == (1, 0, 1, 0) && return -2
    bits == (1, 0, 0, 1) && return -1 - im
    bits == (0, 1, 1, 0) && return -1 - im
    bits == (0, 1, 0, 1) && return -2im
    bits == (0, 0, 1, 1) && return 1 - im
    return 0
end

function gnw_A(bits)
    bits == (1, 1, 0, 0) && return 1 + im
    bits == (1, 0, 1, 0) && return 2
    bits == (1, 0, 0, 1) && return 1 - im
    bits == (0, 1, 1, 0) && return 1 - im
    bits == (0, 1, 0, 1) && return -2im
    bits == (0, 0, 1, 1) && return -1 - im
    return 0
end

# Independent transcription of equations (24)-(30) and the final tensor.
function published_gnw_coefficient(bits, D, g2, z, s)
    i1, j1, i2, j2, j1p, i1p, j2p, i2p = bits
    bar_degree = i1 + i2 + j1p + j2p
    unbar_degree = j1 + j2 + i1p + i2p

    ordering_phase = (
        i1 * (i2p + i1p + j2 + j1) +
        i2 * (i2p + i1p + j2) +
        j1p * (i2p + i1p) +
        j2p * i2p + i1p + i2p
    )
    chemical_power = i2 + i2p - j2 - j2p
    prefactor = (-1)^ordering_phase * s^sum(bits) * z^chemical_power

    local_weight = 0
    if bar_degree == 0 && unbar_degree == 0
        local_weight += D^2 + 2g2
    elseif bar_degree == 1 && unbar_degree == 1
        first_mass_phase = i1p + i2 + i1 + j2
        imaginary_phase = i2p + i2 + j2p + j2
        local_weight -= (-1)^first_mass_phase * im^imaginary_phase * D
        local_weight -= D
    end

    local_weight -= (
        gnw_Abar((i1, i2, j1p, j2p)) *
        gnw_A((j1, j2, i1p, i2p))
    )
    return prefactor * local_weight
end

function gnw_occupation_bits(mask)
    return ntuple(k -> (mask >> (k - 1)) & 1, 8)
end

function gnw_occupied_legs(legs, bits)
    return GrassmannGenerator[
        leg for (leg, occupation) in zip(legs, bits) if occupation == 1
    ]
end

function published_gnw_tensor(legs, D, g2, z, s)
    result = GrassmannExpr(0)
    for mask in 0:255
        bits = gnw_occupation_bits(mask)
        coefficient = published_gnw_coefficient(bits, D, g2, z, s)
        monomial = GrassmannExpr(1)
        for leg in gnw_occupied_legs(legs, bits)
            monomial = monomial * leg
        end
        result = result + coefficient * monomial
    end
    return simplify_coefficients(result)
end

function validate_single_flavor_gnw_tensor()
    model = derive_single_flavor_gnw_tensor()
    D, g2, z, s = model.parameters
    reference = published_gnw_tensor(model.legs, D, g2, z, s)

    @test symbolically_equal(model.tensor, reference)
    @test is_even(model.tensor)
    for mask in 0:255
        bits = gnw_occupation_bits(mask)
        actual = get_coeff(model.tensor, gnw_occupied_legs(model.legs, bits))
        expected = published_gnw_coefficient(bits, D, g2, z, s)
        @test symbolically_equal(GrassmannExpr(actual), GrassmannExpr(expected))
    end
    return model
end

if abspath(PROGRAM_FILE) == @__FILE__
    validate_single_flavor_gnw_tensor()
end
