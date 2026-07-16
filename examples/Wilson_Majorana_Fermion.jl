using GrassmannSymbolics
using Symbolics
using Test

# Wilson-Majorana fermion tensor from
# Ref/Wilson_Majorana_Fermion_and_Two_Flavor_Staggered_GN.pdf, section 2.
#
# The paper separates the initial tensor network into a numerical tensor T and
# a Grassmann tensor G.  This file reconstructs T by direct Berezin integration
# of Eq. (2.9), then packages the coefficients as a GrassmannExpr whose four
# legs are the 4-bit super-indices (i, j, k, l).

function _wm_bit_tuple(mask::Integer, width::Integer=16)
    return ntuple(k -> (mask >> (width - k)) & 1, width)
end

function _wm_chunks(bits)
    length(bits) == 16 || throw(ArgumentError("expected 16 occupation bits"))
    return (Tuple(bits[1:4]), Tuple(bits[5:8]), Tuple(bits[9:12]), Tuple(bits[13:16]))
end

function _wm_monomial(legs, bits)
    result = GrassmannExpr(1)
    for (leg, bit) in zip(legs, bits)
        bit == 0 && continue
        result = result * leg
    end
    return result
end

function _wm_tensor_from_coefficient(legs, coefficient_at_bits)
    result = GrassmannExpr(0)
    for mask in 0:(2^length(legs)-1)
        bits = _wm_bit_tuple(mask, length(legs))
        bit_degree = sum(bits)
        (isodd(bit_degree) || bit_degree > 4) && continue
        result = result + coefficient_at_bits(bits) * _wm_monomial(legs, bits)
    end
    return simplify_coefficients(result)
end

function _wm_selected_legs(legs, bits)
    return GrassmannGenerator[leg for (leg, bit) in zip(legs, bits) if bit == 1]
end

function _wm_generator_powers(factors, bits)
    result = GrassmannExpr(1)
    for (factor, bit) in zip(factors, bits)
        bit == 0 && continue
        result = result * factor
    end
    return result
end

function wilson_majorana_variables()
    @variables m_eta m_chi q

    eta = grassmann(:eta, 1:2)
    chi = grassmann(:chi, 1:2)

    alpha = grassmann(:alpha, 1:4)
    beta = grassmann(:beta, 1:4)
    alphabar = dual.(alpha)
    betabar = dual.(beta)

    eta_tilde = (
        q * (eta[1] + eta[2]),
        q * (-eta[1] + eta[2]),
    )
    chi_tilde = (
        q * (chi[1] + chi[2]),
        q * (-chi[1] + chi[2]),
    )

    physical = (; eta, chi, eta_tilde, chi_tilde)
    legs = (alpha..., beta..., alphabar..., betabar...)
    leg_groups = (collect(alpha), collect(beta), collect(alphabar), collect(betabar))
    parameters = (; m_eta, m_chi, q)
    measure_order = GrassmannGenerator[eta[1], eta[2], chi[1], chi[2]]
    return (; physical, legs, leg_groups, parameters, measure_order)
end

function wilson_majorana_physical_factors(vars, i, j, k, l)
    p = vars.physical
    eta = p.eta
    chi = p.chi
    etat = p.eta_tilde
    chit = p.chi_tilde

    # Eq. (2.9), transcribed in the printed order:
    # eta_2^l4 chi_2^l3 chi_2^l2 eta_2^l1
    # etat_2^k4 chit_2^k3 chit_2^k2 etat_2^k1
    # chi_1^j4 eta_1^j3 chi_1^j2 eta_1^j1
    # chit_1^i4 etat_1^i3 chit_1^i2 etat_1^i1.
    factors = (
        eta[2], chi[2], chi[2], eta[2],
        etat[2], chit[2], chit[2], etat[2],
        chi[1], eta[1], chi[1], eta[1],
        chit[1], etat[1], chit[1], etat[1],
    )
    powers = (l[4], l[3], l[2], l[1], k[4], k[3], k[2], k[1],
              j[4], j[3], j[2], j[1], i[4], i[3], i[2], i[1])
    return _wm_generator_powers(factors, powers)
end

function wilson_majorana_coefficient(vars, bits; simplify::Bool=true)
    i, j, k, l = _wm_chunks(bits)
    p = vars.parameters
    eta = vars.physical.eta
    chi = vars.physical.chi

    bit_degree = sum(bits)
    (isodd(bit_degree) || bit_degree > 4) && return 0

    phase = (-1)^(i[3] + i[4])
    onsite = exp((p.m_eta + 2) * eta[1] * eta[2] +
                 (p.m_chi + 2) * chi[1] * chi[2])
    physical = wilson_majorana_physical_factors(vars, i, j, k, l)
    coefficient = phase * scalar_part(integrate(onsite * physical, vars.measure_order))
    return simplify ? Symbolics.simplify(coefficient; expand=true) : coefficient
end

function derive_wilson_majorana_tensor()
    vars = wilson_majorana_variables()
    tensor = _wm_tensor_from_coefficient(
        vars.legs,
        bits -> wilson_majorana_coefficient(vars, bits),
    )
    return (; tensor, vars..., coefficient=bits -> wilson_majorana_coefficient(vars, bits))
end

function validate_wilson_majorana_tensor()
    @testset "Wilson-Majorana tensor" begin
        model = derive_wilson_majorana_tensor()
        @test is_even(model.tensor)
        @test coefficient_array_shape(model.leg_groups) == (16, 16, 16, 16)
        p = model.parameters
        @test symbolically_equal(
            GrassmannExpr(scalar_part(model.tensor)),
            GrassmannExpr((p.m_eta + 2) * (p.m_chi + 2)),
        )

        for mask in (0, 3, 12, 48, 255, 3840, 61440)
            bits = _wm_bit_tuple(mask)
            occupied = _wm_selected_legs(model.legs, bits)
            actual = get_coeff(model.tensor, occupied)
            expected = model.coefficient(bits)
            @test symbolically_equal(GrassmannExpr(actual), GrassmannExpr(expected))
        end
    end
    return derive_wilson_majorana_tensor()
end

if abspath(PROGRAM_FILE) == @__FILE__
    model = validate_wilson_majorana_tensor()
    println("wilson_majorana_terms_", length(model.tensor.terms))
end
