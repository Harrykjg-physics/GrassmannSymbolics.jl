using GrassmannSymbolics
using Symbolics
using Test

# Two-flavor staggered Gross-Neveu tensor from
# Ref/Wilson_Majorana_Fermion_and_Two_Flavor_Staggered_GN.pdf, section 5.
#
# Eq. (5.2) gives a local 16-bit coefficient tensor.  This file reconstructs it
# by direct Berezin integration over the four staggered components psi_1,...,psi_4.

function _tfsgn_bit_tuple(mask::Integer, width::Integer=16)
    return ntuple(k -> (mask >> (width - k)) & 1, width)
end

function _tfsgn_chunks(bits)
    length(bits) == 16 || throw(ArgumentError("expected 16 occupation bits"))
    return (Tuple(bits[1:4]), Tuple(bits[5:8]), Tuple(bits[9:12]), Tuple(bits[13:16]))
end

function _tfsgn_monomial(legs, bits)
    result = GrassmannExpr(1)
    for (leg, bit) in zip(legs, bits)
        bit == 0 && continue
        result = result * leg
    end
    return result
end

function _tfsgn_tensor_from_coefficient(legs, coefficient_at_bits)
    result = GrassmannExpr(0)
    for mask in 0:(2^length(legs)-1)
        bits = _tfsgn_bit_tuple(mask, length(legs))
        bit_degree = sum(bits)
        (isodd(bit_degree) || bit_degree > 4) && continue
        result = result + coefficient_at_bits(bits) * _tfsgn_monomial(legs, bits)
    end
    return simplify_coefficients(result)
end

function _tfsgn_selected_legs(legs, bits)
    return GrassmannGenerator[leg for (leg, bit) in zip(legs, bits) if bit == 1]
end

function _tfsgn_generator_powers(factors, bits)
    result = GrassmannExpr(1)
    for (factor, bit) in zip(factors, bits)
        bit == 0 && continue
        result = result * factor
    end
    return result
end

function two_flavor_staggered_gn_variables()
    @variables U p_x p_t

    psi = grassmann(:psi, 1:4)
    alpha = grassmann(:sgn_alpha, 1:4)
    beta = grassmann(:sgn_beta, 1:4)
    alphabar = dual.(alpha)
    betabar = dual.(beta)

    physical = (; psi)
    legs = (alpha..., beta..., alphabar..., betabar...)
    leg_groups = (collect(alpha), collect(beta), collect(alphabar), collect(betabar))
    parameters = (; U, p_x, p_t)
    measure_order = GrassmannGenerator[psi[1], psi[2], psi[3], psi[4]]
    return (; physical, legs, leg_groups, parameters, measure_order)
end

function two_flavor_staggered_gn_physical_factors(vars, i, j, k, l)
    psi = vars.physical.psi

    # Eq. (5.2), transcribed in the printed order:
    # psi_4^l4 psi_3^l3 psi_2^l2 psi_1^l1
    # psi_4^k4 ... psi_1^k1
    # psi_4^j4 ... psi_1^j1
    # psi_4^i4 ... psi_1^i1.
    factors = (
        psi[4], psi[3], psi[2], psi[1],
        psi[4], psi[3], psi[2], psi[1],
        psi[4], psi[3], psi[2], psi[1],
        psi[4], psi[3], psi[2], psi[1],
    )
    powers = (l[4], l[3], l[2], l[1], k[4], k[3], k[2], k[1],
              j[4], j[3], j[2], j[1], i[4], i[3], i[2], i[1])
    return _tfsgn_generator_powers(factors, powers)
end

function two_flavor_staggered_gn_coefficient(vars, bits; simplify::Bool=true)
    i, j, k, l = _tfsgn_chunks(bits)
    psi = vars.physical.psi
    p = vars.parameters

    bit_degree = sum(bits)
    (isodd(bit_degree) || bit_degree > 4) && return 0

    # Repeated physical psi fields vanish before the onsite four-fermi factor
    # can contribute.
    selected = Int[]
    append!(selected, [5 - a for a in 1:4 if l[a] == 1])
    append!(selected, [5 - a for a in 1:4 if k[a] == 1])
    append!(selected, [5 - a for a in 1:4 if j[a] == 1])
    append!(selected, [5 - a for a in 1:4 if i[a] == 1])
    length(unique(selected)) == length(selected) || return 0

    staggered_weight = (p.p_x * (1 // 2))^sum(i) * (p.p_t * (1 // 2))^sum(j)
    onsite = exp(p.U * psi[1] * psi[2] * psi[3] * psi[4])
    physical = two_flavor_staggered_gn_physical_factors(vars, i, j, k, l)
    coefficient = staggered_weight *
        scalar_part(integrate(onsite * physical, vars.measure_order))
    return simplify ? Symbolics.simplify(coefficient; expand=true) : coefficient
end

function derive_two_flavor_staggered_gn_tensor()
    vars = two_flavor_staggered_gn_variables()
    tensor = _tfsgn_tensor_from_coefficient(
        vars.legs,
        bits -> two_flavor_staggered_gn_coefficient(vars, bits),
    )
    return (; tensor, vars..., coefficient=bits -> two_flavor_staggered_gn_coefficient(vars, bits))
end

function validate_two_flavor_staggered_gn_tensor()
    @testset "Two-flavor staggered Gross-Neveu tensor" begin
        model = derive_two_flavor_staggered_gn_tensor()
        @test is_even(model.tensor)
        @test coefficient_array_shape(model.leg_groups) == (16, 16, 16, 16)
        @test symbolically_equal(
            GrassmannExpr(scalar_part(model.tensor)),
            GrassmannExpr(model.parameters.U),
        )
        @test symbolically_equal(
            GrassmannExpr(model.coefficient((1, 1, 1, 1, ntuple(_ -> 0, 12)...))),
            GrassmannExpr((model.parameters.p_x * (1 // 2))^4),
        )

        for mask in (0, 1, 15, 240, 3840, 61440, 65535)
            bits = _tfsgn_bit_tuple(mask)
            occupied = _tfsgn_selected_legs(model.legs, bits)
            actual = get_coeff(model.tensor, occupied)
            expected = model.coefficient(bits)
            @test symbolically_equal(GrassmannExpr(actual), GrassmannExpr(expected))
        end
    end
    return derive_two_flavor_staggered_gn_tensor()
end

if abspath(PROGRAM_FILE) == @__FILE__
    model = validate_two_flavor_staggered_gn_tensor()
    println("two_flavor_staggered_gn_terms_", length(model.tensor.terms))
end
