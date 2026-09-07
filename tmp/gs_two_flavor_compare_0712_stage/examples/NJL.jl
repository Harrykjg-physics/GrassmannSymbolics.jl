using GrassmannSymbolics
using Symbolics
using Test

# Four-dimensional lattice NJL model of JHEP 01 (2021) 121.
# This file implements equations (2.11)--(2.15).  Each direction carries
# a 3-bit index: channel 1 is the forward hopping split, channel 2 is the
# backward hopping split, and channel 3 is the four-fermi split.  The order of
# directions follows the paper: t, x, y, z, t0, x0, y0, z0.

struct NJLBits
    t::NTuple{3,Int}
    x::NTuple{3,Int}
    y::NTuple{3,Int}
    z::NTuple{3,Int}
    tp::NTuple{3,Int}
    xp::NTuple{3,Int}
    yp::NTuple{3,Int}
    zp::NTuple{3,Int}
end

function njl_bits(values)
    length(values) == 24 || throw(ArgumentError(string(:NJL_bits_need_24_entries)))
    return NJLBits(
        Tuple(values[1:3]), Tuple(values[4:6]), Tuple(values[7:9]),
        Tuple(values[10:12]), Tuple(values[13:15]), Tuple(values[16:18]),
        Tuple(values[19:21]), Tuple(values[22:24]),
    )
end

function njl_bits_from_mask(mask)
    return njl_bits(ntuple(k -> (mask >> (k - 1)) & 1, 24))
end

function njl_all_direction_triples(bits::NJLBits)
    return (bits.t, bits.x, bits.y, bits.z, bits.tp, bits.xp, bits.yp, bits.zp)
end

function njl_unprimed_triples(bits::NJLBits)
    return (bits.t, bits.x, bits.y, bits.z)
end

function njl_primed_triples(bits::NJLBits)
    return (bits.tp, bits.xp, bits.yp, bits.zp)
end

δ(a, b) = a == b ? 1 : 0

function njl_bar_degree(bits::NJLBits)
    triples = njl_all_direction_triples(bits)
    return sum(triple[3] for triple in triples) +
        sum(triple[1] for triple in njl_unprimed_triples(bits)) +
        sum(triple[2] for triple in njl_primed_triples(bits))
end

function njl_chi_degree(bits::NJLBits)
    triples = njl_all_direction_triples(bits)
    return sum(triple[3] for triple in triples) +
        sum(triple[2] for triple in njl_unprimed_triples(bits)) +
        sum(triple[1] for triple in njl_primed_triples(bits))
end

njl_bar_delta(bits::NJLBits, q) = δ(njl_bar_degree(bits), q)
njl_chi_delta(bits::NJLBits, q) = δ(njl_chi_degree(bits), q)

function njl_staggered_phase(bits::NJLBits, n1, n2, n3)
    t1, t2 = bits.t[1], bits.t[2]
    y1, y2 = bits.y[1], bits.y[2]
    z1, z2 = bits.z[1], bits.z[2]
    return n1 * (y1 + y2 + z1 + z2 + t1 + t2) +
        n2 * (z1 + z2 + t1 + t2) + n3 * (t1 + t2)
end

function njl_hopping_degree(bits::NJLBits)
    return sum(triple[1] + triple[2] for triple in njl_all_direction_triples(bits))
end

function njl_fourfermi_degree(bits::NJLBits)
    return sum(triple[3] for triple in njl_all_direction_triples(bits))
end

function njl_temporal_chemical_degree(bits::NJLBits)
    return bits.t[1] - bits.t[2] + bits.tp[1] - bits.tp[2]
end

function njl_I_coefficient(bits::NJLBits, m, r, u, n1, n2, n3)
    local_integral = -m * njl_bar_delta(bits, 0) * njl_chi_delta(bits, 0) +
        njl_bar_delta(bits, 1) * njl_chi_delta(bits, 1)
    return (-1)^njl_staggered_phase(bits, n1, n2, n3) *
        (-1 // 2)^njl_hopping_degree(bits) *
        r^njl_fourfermi_degree(bits) *
        u^njl_temporal_chemical_degree(bits) *
        local_integral
end

function njl_S_phase(bits::NJLBits)
    t1, t2 = bits.t[1], bits.t[2]
    x1, x2 = bits.x[1], bits.x[2]
    y1, y2 = bits.y[1], bits.y[2]
    z1, z2 = bits.z[1], bits.z[2]
    tp1, tp2 = bits.tp[1], bits.tp[2]
    xp1, xp2 = bits.xp[1], bits.xp[2]
    yp1, yp2 = bits.yp[1], bits.yp[2]
    zp1, zp2 = bits.zp[1], bits.zp[2]

    current_phase = t1 * (t2 + x2 + y2 + z2) +
        x1 * (x2 + y2 + z2) + y1 * (y2 + z2) + z1 * z2
    primed_phase = tp2 * (tp1 + xp1 + yp1 + zp1) +
        xp2 * (xp1 + yp1 + zp1) + yp2 * (yp1 + zp1) + zp2 * zp1
    bridge_phase = (t1 + t2 + x1 + x2 + y1 + y2 + z1 + z2) *
        (tp1 + xp1 + yp1 + zp1)
    return current_phase + primed_phase + bridge_phase
end

njl_S_sign(bits::NJLBits) = (-1)^njl_S_phase(bits)

function njl_variables()
    @variables m r u
    φt = grassmann(:φt)
    ψt = grassmann(:ψt)
    φx = grassmann(:φx)
    ψx = grassmann(:ψx)
    φy = grassmann(:φy)
    ψy = grassmann(:ψy)
    φz = grassmann(:φz)
    ψz = grassmann(:ψz)
    φbart = dual(φt)
    ψbart = dual(ψt)
    φbarx = dual(φx)
    ψbarx = dual(ψx)
    φbary = dual(φy)
    ψbary = dual(ψy)
    φbarz = dual(φz)
    ψbarz = dual(ψz)
    legs = (; φt, ψt, φx, ψx, φy, ψy, φz, ψz,
        ψbart, φbart, ψbarx, φbarx, ψbary, φbary, ψbarz, φbarz)
    parameters = (; m, r, u)
    return (; legs, parameters)
end

function _multiply_if_power(result, factor, power)
    power == 0 && return result
    return result * factor
end

function njl_G_monomial(bits::NJLBits, legs)
    result = GrassmannExpr(1)
    bilinears = (
        legs.φbart * legs.φt, legs.ψbart * legs.ψt,
        legs.φbarx * legs.φx, legs.ψbarx * legs.ψx,
        legs.φbary * legs.φy, legs.ψbary * legs.ψy,
        legs.φbarz * legs.φz, legs.ψbarz * legs.ψz,
    )
    powers = (bits.t[1], bits.t[2], bits.x[1], bits.x[2],
        bits.y[1], bits.y[2], bits.z[1], bits.z[2])
    for (factor, power) in zip(bilinears, powers)
        result = _multiply_if_power(result, factor, power)
    end
    return result
end

function njl_raw_G_monomial(bits::NJLBits, legs)
    result = GrassmannExpr(1)
    bilinears = (
        legs.φbarx * legs.φx, legs.ψbarx * legs.ψx,
        legs.φbary * legs.φy, legs.ψbary * legs.ψy,
        legs.φbarz * legs.φz, legs.ψbarz * legs.ψz,
        legs.φbart * legs.φt, legs.ψbart * legs.ψt,
    )
    powers = (bits.x[1], bits.x[2], bits.y[1], bits.y[2],
        bits.z[1], bits.z[2], bits.t[1], bits.t[2])
    for (factor, power) in zip(bilinears, powers)
        result = _multiply_if_power(result, factor, power)
    end
    return result
end

function njl_G_measure_order(bits::NJLBits, legs)
    powers_and_legs = (
        (bits.t[1], legs.φt), (bits.t[2], legs.ψt),
        (bits.x[1], legs.φx), (bits.x[2], legs.ψx),
        (bits.y[1], legs.φy), (bits.y[2], legs.ψy),
        (bits.z[1], legs.φz), (bits.z[2], legs.ψz),
        (bits.tp[2], legs.ψbart), (bits.tp[1], legs.φbart),
        (bits.xp[2], legs.ψbarx), (bits.xp[1], legs.φbarx),
        (bits.yp[2], legs.ψbary), (bits.yp[1], legs.φbary),
        (bits.zp[2], legs.ψbarz), (bits.zp[1], legs.φbarz),
    )
    order = GrassmannGenerator[]
    for (power, leg) in powers_and_legs
        power == 1 && push!(order, leg)
    end
    return order
end

function njl_T_entry(bits::NJLBits, model; n1=0, n2=0, n3=0)
    p = model.parameters
    coeff = njl_I_coefficient(bits, p.m, p.r, p.u, n1, n2, n3) *
        njl_S_sign(bits)
    return simplify_coefficients(coeff * njl_G_monomial(bits, model.legs))
end

function njl_projected_entry_without_measure_sign(bits::NJLBits, model; n1=0, n2=0, n3=0)
    p = model.parameters
    coeff = njl_I_coefficient(bits, p.m, p.r, p.u, n1, n2, n3)
    return simplify_coefficients(coeff * njl_G_monomial(bits, model.legs))
end

function njl_direct_local_integral(bits::NJLBits, m, r, u, n1, n2, n3)
    χ = grassmann(:χ)
    χbar = dual(χ)
    η = (1, (-1)^n1, (-1)^(n1 + n2), (-1)^(n1 + n2 + n3))
    temporal = (1, 1, 1, u)
    inverse_temporal = (1, 1, 1, u^-1)
    current = (bits.x, bits.y, bits.z, bits.t)
    previous = (bits.xp, bits.yp, bits.zp, bits.tp)

    integrand = exp(-m * χbar * χ)
    for direction in 1:4
        local_now = current[direction]
        local_prev = previous[direction]
        integrand = integrand * (((-η[direction] * temporal[direction] / 2) * χbar) ^ local_now[1])
        integrand = integrand * (((-η[direction] * inverse_temporal[direction] / 2) * χ) ^ local_now[2])
        integrand = integrand * ((r * χbar * χ) ^ local_now[3])
        integrand = integrand * (((-temporal[direction] / 2) * χ) ^ local_prev[1])
        integrand = integrand * (((-inverse_temporal[direction] / 2) * χbar) ^ local_prev[2])
        integrand = integrand * ((r * χbar * χ) ^ local_prev[3])
    end
    return scalar_part(integrate(integrand, [χ, χbar]))
end

function njl_normal_ordered_local_integral(bits::NJLBits, m, r, u, n1, n2, n3)
    χ = grassmann(:χ)
    χbar = dual(χ)
    prefactor = (-1)^njl_staggered_phase(bits, n1, n2, n3) *
        (-1 // 2)^njl_hopping_degree(bits) *
        r^njl_fourfermi_degree(bits) *
        u^njl_temporal_chemical_degree(bits)
    integrand = exp(-m * χbar * χ) *
        (χbar ^ njl_bar_degree(bits)) * (χ ^ njl_chi_degree(bits))
    return prefactor * scalar_part(integrate(integrand, [χ, χbar]))
end

function njl_nonzero_candidate_bits()
    candidates = NJLBits[]
    for mask in 0:(2^24 - 1)
        bits = njl_bits_from_mask(mask)
        bar_degree = njl_bar_degree(bits)
        chi_degree = njl_chi_degree(bits)
        if ((bar_degree == 0) && (chi_degree == 0)) ||
                ((bar_degree == 1) && (chi_degree == 1))
            push!(candidates, bits)
        end
    end
    return candidates
end

function validate_njl_tensor()
    model = njl_variables()
    p = model.parameters
    candidates = njl_nonzero_candidate_bits()
    @test length(candidates) == 73

    for n1 in 0:1, n2 in 0:1, n3 in 0:1
        for bits in candidates
            direct = njl_normal_ordered_local_integral(bits, p.m, p.r, p.u, n1, n2, n3)
            published = njl_I_coefficient(bits, p.m, p.r, p.u, n1, n2, n3)
            @test symbolically_equal(GrassmannExpr(direct), GrassmannExpr(published))
            @test is_even(njl_G_monomial(bits, model.legs))
        end
    end

    @test njl_S_sign(njl_bits((1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))) == -1
    @test njl_S_sign(njl_bits((0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))) == 1

    zero_samples = (1, 7, 63, 511, 4095, 32767, 262143, 2097151, 16777215)
    for mask in zero_samples
        bits = njl_bits_from_mask(mask)
        if !(bits in candidates)
        @test njl_I_coefficient(bits, p.m, p.r, p.u, 1, 1, 1) == 0
        end
    end

    empty_bits = njl_bits(ntuple(_ -> 0, 24))
    @test symbolically_equal(
        GrassmannExpr(njl_I_coefficient(empty_bits, p.m, p.r, p.u, 0, 0, 0)),
        GrassmannExpr(-p.m),
    )
    @test njl_S_sign(empty_bits) == 1
    @test isempty(njl_G_measure_order(empty_bits, model.legs))
    return (; model, candidates)
end

if abspath(PROGRAM_FILE) == @__FILE__
    validate_njl_tensor()
end
