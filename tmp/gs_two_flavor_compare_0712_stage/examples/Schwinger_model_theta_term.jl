using GrassmannSymbolics
using Symbolics
using Test

# Fermionic tensor in appendix C. Here uν=exp(i*pi*aν), p=(-1)^n1,
# and s=1/sqrt(2). The reference uses the measure d(chibar)d(chi).
function derive_schwinger_fermion_tensor()
    @variables m u1 u2 s p

    χ = grassmann(:χ)
    χbar = dual(χ)
    ζ1, ξ1, ζ2, ξ2 = (
        grassmann(:ζ1), grassmann(:ξ1), grassmann(:ζ2), grassmann(:ξ2)
    )
    ζbar1, ξbar1, ζbar2, ξbar2 = dual.((ζ1, ξ1, ζ2, ξ2))

    bond_exponents = [
        s * u1 * χbar * ζ1,
        s * χ * ζbar1,
        (s / u1) * χ * ξ1,
        s * χbar * ξbar1,
        s * p * u2 * χbar * ζ2,
        s * χ * ζbar2,
        (s * p / u2) * χ * ξ2,
        s * χbar * ξbar2,
    ]
    tensor = local_grassmann_tensor(
        m * χbar * χ,
        bond_exponents,
        [χbar, χ],
    )

    legs = (ζ1, ξ1, ζ2, ξ2, ξbar1, ζbar1, ξbar2, ζbar2)
    parameters = (; m, u1, u2, s, p)
    return (; tensor, legs, parameters)
end

# Bosonic plaquette factor in equation (3.6), before index contractions.
function schwinger_gauge_plaquette_factor(flux, β, θ, quadrature_weights)
    quadrature_factor = sqrt(prod(quadrature_weights)) / 4
    action_factor = exp(
        β * cos(pi * flux) +
        (θ / (2pi)) * log(exp(im * pi * flux))
    )
    return quadrature_factor * action_factor
end

function published_schwinger_fermion_coefficient(bits, m, u1, u2, s, p)
    i1, j1, i2, j2, j1p, i1p, j2p, i2p = bits
    bar_degree = i1 + i2 + j1p + j2p
    unbar_degree = i1p + i2p + j1 + j2

    local_weight = if bar_degree == 1 && unbar_degree == 1
        1
    elseif bar_degree == 0 && unbar_degree == 0
        m
    else
        0
    end
    phase = (
        j1 * (i2 + j1p + j2p) +
        j2 * (j1p + j2p) +
        i1p * j2p
    )
    gauge_factor = u1^(i1 - j1) * u2^(i2 - j2)
    staggered_factor = p^(i2 + j2)
    return (
        local_weight * (-1)^phase * s^sum(bits) *
        gauge_factor * staggered_factor
    )
end

function schwinger_occupation_bits(mask)
    return ntuple(k -> (mask >> (k - 1)) & 1, 8)
end

function schwinger_occupied_legs(legs, bits)
    return GrassmannGenerator[
        leg for (leg, occupation) in zip(legs, bits) if occupation == 1
    ]
end

function published_schwinger_fermion_tensor(model)
    m, u1, u2, s, p = model.parameters
    result = GrassmannExpr(0)
    for mask in 0:255
        bits = schwinger_occupation_bits(mask)
        coefficient = published_schwinger_fermion_coefficient(
            bits, m, u1, u2, s, p,
        )
        monomial = GrassmannExpr(1)
        for leg in schwinger_occupied_legs(model.legs, bits)
            monomial = monomial * leg
        end
        result = result + coefficient * monomial
    end
    return simplify_coefficients(result)
end

function validate_schwinger_fermion_tensor()
    model = derive_schwinger_fermion_tensor()
    reference = published_schwinger_fermion_tensor(model)
    m, u1, u2, s, p = model.parameters

    @test symbolically_equal(model.tensor, reference)
    @test is_even(model.tensor)
    for mask in 0:255
        bits = schwinger_occupation_bits(mask)
        occupied = schwinger_occupied_legs(model.legs, bits)
        actual = get_coeff(model.tensor, occupied)
        expected = published_schwinger_fermion_coefficient(
            bits, m, u1, u2, s, p,
        )
        @test symbolically_equal(GrassmannExpr(actual), GrassmannExpr(expected))
    end
    return model
end

if abspath(PROGRAM_FILE) == @__FILE__
    validate_schwinger_fermion_tensor()
end
