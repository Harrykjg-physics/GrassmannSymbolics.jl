using GrassmannSymbolics
using Symbolics
using Test

# Appendix A of the (1+1)-dimensional Hubbard TRG paper.
# Parameters: q = sqrt(t * epsilon), h = mu * epsilon + 1, Ue = U * epsilon.
# Leg order follows equation (4):
#   Ψσ=(ησ↑,ησ↓,ζσ↑,ζσ↓), Ψτ=(ητ↑,ητ↓),
#   Ψbarτ=(ηbarτ↓,ηbarτ↑), Ψbarσ=(ζbarσ↓,ζbarσ↑,ηbarσ↓,ηbarσ↑).

struct HubbardBits
    iσu::Int
    iσd::Int
    jσu::Int
    jσd::Int
    iτu::Int
    iτd::Int
    iσpu::Int
    iσpd::Int
    jσpu::Int
    jσpd::Int
    iτpu::Int
    iτpd::Int
end

function hubbard_bits_from_leg_bits(bits)
    length(bits) == 12 || throw(ArgumentError(string(:Hubbard_needs_12_bits)))
    return HubbardBits(
        bits[1], bits[2], bits[3], bits[4], bits[5], bits[6],
        bits[12], bits[11], bits[10], bits[9], bits[8], bits[7],
    )
end

function hubbard_leg_bits(mask)
    return ntuple(k -> (mask >> (k - 1)) & 1, 12)
end

δ(a, b) = a == b ? 1 : 0

function hubbard_variables()
    @variables q h Ue

    ψu = grassmann(:ψ, index=1)
    ψd = grassmann(:ψ, index=2)
    ψbaru = dual(ψu)
    ψbard = dual(ψd)

    ησu = grassmann(:ησ, index=1)
    ησd = grassmann(:ησ, index=2)
    ζσu = grassmann(:ζσ, index=1)
    ζσd = grassmann(:ζσ, index=2)
    ητu = grassmann(:ητ, index=1)
    ητd = grassmann(:ητ, index=2)

    ηbarσu = dual(ησu)
    ηbarσd = dual(ησd)
    ζbarσu = dual(ζσu)
    ζbarσd = dual(ζσd)
    ηbarτu = dual(ητu)
    ηbarτd = dual(ητd)

    original = (; ψu, ψd, ψbaru, ψbard)
    legs = (
        ησu, ησd, ζσu, ζσd, ητu, ητd,
        ηbarτd, ηbarτu, ζbarσd, ζbarσu, ηbarσd, ηbarσu,
    )
    named_legs = (;
        ησu, ησd, ζσu, ζσd, ητu, ητd,
        ηbarτd, ηbarτu, ζbarσd, ζbarσu, ηbarσd, ηbarσu,
    )
    parameters = (; q, h, Ue)
    return (; original, legs, named_legs, parameters)
end

function hubbard_spatial_degree(bits::HubbardBits)
    return bits.iσu + bits.iσd + bits.jσu + bits.jσd +
        bits.iσpu + bits.iσpd + bits.jσpu + bits.jσpd
end

function hubbard_forward_spatial_degree(bits::HubbardBits)
    return bits.iσu + bits.iσd
end

function hubbard_delta_arguments(bits::HubbardBits)
    au = bits.iτu + bits.iσu + bits.jσpu
    bu = bits.iτpu + bits.iσpu + bits.jσu
    ad = bits.iτd + bits.iσd + bits.jσpd
    bd = bits.iτpd + bits.iσpd + bits.jσd
    return (; au, bu, ad, bd)
end

function _permutation_sign_between_orders(raw, target)
    filtered_raw = [generator for generator in raw if generator !== nothing]
    filtered_target = [generator for generator in target if generator !== nothing]
    length(filtered_raw) == length(filtered_target) || return 0
    sort(filtered_raw) == sort(filtered_target) || return 0
    positions = [findfirst(==(generator), filtered_target) for generator in filtered_raw]
    inversions = 0
    for i in 1:length(positions), j in (i + 1):length(positions)
        positions[i] > positions[j] && (inversions += 1)
    end
    return iseven(inversions) ? 1 : -1
end

function hubbard_R_sign(bits::HubbardBits, named_legs)
    raw = Any[
        bits.iτu == 1 ? named_legs.ητu : nothing,
        bits.iσu == 1 ? named_legs.ησu : nothing,
        bits.jσu == 1 ? named_legs.ζσu : nothing,
        bits.iτpu == 1 ? named_legs.ηbarτu : nothing,
        bits.iσpu == 1 ? named_legs.ηbarσu : nothing,
        bits.jσpu == 1 ? named_legs.ζbarσu : nothing,
        bits.iτd == 1 ? named_legs.ητd : nothing,
        bits.iσd == 1 ? named_legs.ησd : nothing,
        bits.jσd == 1 ? named_legs.ζσd : nothing,
        bits.iτpd == 1 ? named_legs.ηbarτd : nothing,
        bits.iσpd == 1 ? named_legs.ηbarσd : nothing,
        bits.jσpd == 1 ? named_legs.ζbarσd : nothing,
    ]
    target = Any[
        bits.iσu == 1 ? named_legs.ησu : nothing,
        bits.iσd == 1 ? named_legs.ησd : nothing,
        bits.jσu == 1 ? named_legs.ζσu : nothing,
        bits.jσd == 1 ? named_legs.ζσd : nothing,
        bits.iτu == 1 ? named_legs.ητu : nothing,
        bits.iτd == 1 ? named_legs.ητd : nothing,
        bits.iτpd == 1 ? named_legs.ηbarτd : nothing,
        bits.iτpu == 1 ? named_legs.ηbarτu : nothing,
        bits.jσpd == 1 ? named_legs.ζbarσd : nothing,
        bits.jσpu == 1 ? named_legs.ζbarσu : nothing,
        bits.iσpd == 1 ? named_legs.ηbarσd : nothing,
        bits.iσpu == 1 ? named_legs.ηbarσu : nothing,
    ]
    return _permutation_sign_between_orders(raw, target)
end

function published_hubbard_coefficient(bits::HubbardBits, q, h, Ue, named_legs)
    ψu = grassmann(:ψ, index=1)
    ψd = grassmann(:ψ, index=2)
    ψbaru = dual(ψu)
    ψbard = dual(ψd)
    nup = ψbaru * ψu
    ndn = ψbard * ψd
    integrand = exp(-(Ue * nup * ndn - h * (nup + ndn)))

    factors = (
        q * ψbaru * named_legs.ησu,
        q * ψbard * named_legs.ησd,
        q * named_legs.ζσu * ψu,
        q * named_legs.ζσd * ψd,
        -ψbaru * named_legs.ητu,
        -ψbard * named_legs.ητd,
        named_legs.ηbarτd * ψd,
        named_legs.ηbarτu * ψu,
        -q * ψbard * named_legs.ζbarσd,
        -q * ψbaru * named_legs.ζbarσu,
        q * named_legs.ηbarσd * ψd,
        q * named_legs.ηbarσu * ψu,
    )
    powers = (
        bits.iσu, bits.iσd, bits.jσu, bits.jσd, bits.iτu, bits.iτd,
        bits.iτpd, bits.iτpu, bits.jσpd, bits.jσpu, bits.iσpd, bits.iσpu,
    )
    for (factor, power) in zip(factors, powers)
        power == 1 && (integrand = integrand * factor)
    end
    projected = integrate(integrand, [ψbaru, ψu, ψbard, ψd])
    occupied = GrassmannGenerator[]
    leg_bits = (
        bits.iσu, bits.iσd, bits.jσu, bits.jσd, bits.iτu, bits.iτd,
        bits.iτpd, bits.iτpu, bits.jσpd, bits.jσpu, bits.iσpd, bits.iσpu,
    )
    leg_order = (
        named_legs.ησu, named_legs.ησd, named_legs.ζσu, named_legs.ζσd,
        named_legs.ητu, named_legs.ητd, named_legs.ηbarτd, named_legs.ηbarτu,
        named_legs.ζbarσd, named_legs.ζbarσu, named_legs.ηbarσd, named_legs.ηbarσu,
    )
    for (leg, power) in zip(leg_order, leg_bits)
        power == 1 && push!(occupied, leg)
    end
    return simplify_coefficients(GrassmannExpr(get_coeff(projected, occupied))) |> scalar_part
end

function derive_hubbard_tensor()
    model = hubbard_variables()
    p = model.parameters
    o = model.original
    l = model.named_legs

    nup = o.ψbaru * o.ψu
    ndn = o.ψbard * o.ψd
    onsite_action = p.Ue * nup * ndn - p.h * (nup + ndn)

    bond_exponents = [
        p.q * o.ψbaru * l.ησu,
        p.q * o.ψbard * l.ησd,
        p.q * l.ζσu * o.ψu,
        p.q * l.ζσd * o.ψd,
        -o.ψbaru * l.ητu,
        -o.ψbard * l.ητd,
        l.ηbarτu * o.ψu,
        l.ηbarτd * o.ψd,
        -p.q * o.ψbaru * l.ζbarσu,
        -p.q * o.ψbard * l.ζbarσd,
        p.q * l.ηbarσu * o.ψu,
        p.q * l.ηbarσd * o.ψd,
    ]
    tensor = local_grassmann_tensor(
        onsite_action,
        bond_exponents,
        [o.ψbaru, o.ψu, o.ψbard, o.ψd],
    )
    return merge(model, (; tensor))
end

function published_hubbard_tensor(model)
    p = model.parameters
    result = GrassmannExpr(0)
    for mask in 0:4095
        leg_bits = hubbard_leg_bits(mask)
        bits = hubbard_bits_from_leg_bits(leg_bits)
        coefficient = published_hubbard_coefficient(bits, p.q, p.h, p.Ue, model.named_legs)
        monomial = GrassmannExpr(1)
        for (leg, occupation) in zip(model.legs, leg_bits)
            occupation == 1 && (monomial = monomial * leg)
        end
        result = result + coefficient * monomial
    end
    return simplify_coefficients(result)
end

function validate_hubbard_tensor()
    model = derive_hubbard_tensor()
    reference = published_hubbard_tensor(model)
    p = model.parameters

    @test symbolically_equal(model.tensor, reference)
    @test is_even(model.tensor)
    @test length(model.tensor.terms) == length(reference.terms)

    for mask in 0:4095
        leg_bits = hubbard_leg_bits(mask)
        bits = hubbard_bits_from_leg_bits(leg_bits)
        occupied = GrassmannGenerator[
            leg for (leg, occupation) in zip(model.legs, leg_bits) if occupation == 1
        ]
        actual = get_coeff(model.tensor, occupied)
        expected = published_hubbard_coefficient(bits, p.q, p.h, p.Ue, model.named_legs)
        @test symbolically_equal(GrassmannExpr(actual), GrassmannExpr(expected))
    end
    return model
end

if abspath(PROGRAM_FILE) == @__FILE__
    validate_hubbard_tensor()
end
