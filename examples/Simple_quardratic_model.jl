using GrassmannSymbolics
using Symbolics
using Test

# Equation (14) of J. Phys.: Condens. Matter 36 (2024) 343002.
# The parameter q is sqrt(t), where t is the hopping coefficient.
function derive_simple_quadratic_tensor()
    @variables m q

    ψ = grassmann(:ψ)
    ψbar = dual(ψ)
    Φ1, Ψ1, Φ2, Ψ2 = (
        grassmann(:Φ1), grassmann(:Ψ1), grassmann(:Φ2), grassmann(:Ψ2)
    )
    Φbar1, Ψbar1, Φbar2, Ψbar2 = dual.((Φ1, Ψ1, Φ2, Ψ2))

    bond_exponents = [
        q * ψbar * Φbar1,
        q * ψbar * Ψ1,
        q * ψ * Φ1,
        -q * ψ * Ψbar1,
        q * ψbar * Φbar2,
        q * ψbar * Ψ2,
        q * ψ * Φ2,
        -q * ψ * Ψbar2,
    ]
    tensor = local_grassmann_tensor(
        m * ψbar * ψ,
        bond_exponents,
        [ψ, ψbar],
    )

    # This is the precise leg order in equations (15)-(16).
    legs = (Φ1, Ψ1, Φ2, Ψ2, Ψbar1, Φbar1, Ψbar2, Φbar2)
    parameters = (; m, q)
    return (; tensor, legs, parameters)
end

# Independent transcription of the coefficient tensor in equation (16).
function published_simple_coefficient(bits, m, q)
    i1, j1, i2, j2, j1p, i1p, j2p, i2p = bits
    left_degree = i1 + i2 + j1p + j2p
    right_degree = i1p + i2p + j1 + j2

    local_weight = if left_degree == 0 && right_degree == 0
        -m
    elseif left_degree == 1 && right_degree == 1
        1
    else
        0
    end

    phase = (
        j1p + j2p +
        j1 * (i2 + j1p + j2p) +
        j2 * (j1p + j2p) +
        i1p * j2p
    )
    return (-1)^phase * q^sum(bits) * local_weight
end

function occupation_bits(mask, number_of_legs)
    return ntuple(k -> (mask >> (k - 1)) & 1, number_of_legs)
end

function occupied_legs(legs, bits)
    return GrassmannGenerator[
        leg for (leg, occupation) in zip(legs, bits) if occupation == 1
    ]
end

function published_simple_tensor(legs, m, q)
    result = GrassmannExpr(0)
    for mask in 0:(2^length(legs) - 1)
        bits = occupation_bits(mask, length(legs))
        coefficient = published_simple_coefficient(bits, m, q)
        monomial = GrassmannExpr(1)
        for leg in occupied_legs(legs, bits)
            monomial = monomial * leg
        end
        result = result + coefficient * monomial
    end
    return simplify_coefficients(result)
end

function contract_periodic_simple_one_site(model)
    Φ1, Ψ1, Φ2, Ψ2, Ψbar1, Φbar1, Ψbar2, Φbar2 = model.legs
    pairs = (
        (Φbar1, Φ1),
        (Ψbar1, Ψ1),
        (Φbar2, Φ2),
        (Ψbar2, Ψ2),
    )
    return scalar_part(contract_grassmann((model.tensor,), pairs))
end

function validate_simple_quadratic_tensor()
    model = derive_simple_quadratic_tensor()
    m, q = model.parameters
    reference = published_simple_tensor(model.legs, m, q)

    @test symbolically_equal(model.tensor, reference)
    @test length(model.tensor.terms) == 17
    @test is_even(model.tensor)

    for mask in 0:(2^length(model.legs) - 1)
        bits = occupation_bits(mask, length(model.legs))
        actual = get_coeff(model.tensor, occupied_legs(model.legs, bits))
        expected = published_simple_coefficient(bits, m, q)
        @test symbolically_equal(GrassmannExpr(actual), GrassmannExpr(expected))
    end

    closure = contract_periodic_simple_one_site(model)
    @test symbolically_equal(
        GrassmannExpr(closure),
        GrassmannExpr(4 * q^2 - m),
    )
    return model
end

if abspath(PROGRAM_FILE) == @__FILE__
    validate_simple_quadratic_tensor()
end
