using GrassmannSymbolics
using Symbolics
using Test

# Local tensor entries for Eq. (1) of
# "Phase diagram of a lattice fermion model with symmetric mass generation",
# specialized to UB = 0.
#
# Each site has six incident lattice sides, ordered as (t, x, y, tp, xp, yp).
# Without the UB bond four-fermi term, each side carries only four hopping
# occupation bits:
#
#   (u_forward, u_backward, d_forward, d_backward).
#
# The onsite UI term is retained as exp(UI * ubar*u*dbar*d) in exp(-S).

struct SMGNoUBSideBits
    u_forward::Int
    u_backward::Int
    d_forward::Int
    d_backward::Int
end

struct SMGNoUBBits
    t::SMGNoUBSideBits
    x::SMGNoUBSideBits
    y::SMGNoUBSideBits
    tp::SMGNoUBSideBits
    xp::SMGNoUBSideBits
    yp::SMGNoUBSideBits
end

function _smg_no_ub_check_bit(bit)
    bit in (0, 1) ||
        throw(ArgumentError("SMG no-UB occupation bits must be 0 or 1"))
    return Int(bit)
end

function smg_no_ub_side_bits(values)
    length(values) == 4 ||
        throw(ArgumentError("each SMG no-UB side needs 4 occupation bits"))
    checked = Tuple(_smg_no_ub_check_bit(bit) for bit in values)
    return SMGNoUBSideBits(checked...)
end

function smg_no_ub_bits(values)
    length(values) == 24 ||
        throw(ArgumentError("SMG no-UB local tensor entries need 24 occupation bits"))
    return SMGNoUBBits(
        smg_no_ub_side_bits(values[1:4]),
        smg_no_ub_side_bits(values[5:8]),
        smg_no_ub_side_bits(values[9:12]),
        smg_no_ub_side_bits(values[13:16]),
        smg_no_ub_side_bits(values[17:20]),
        smg_no_ub_side_bits(values[21:24]),
    )
end

smg_no_ub_bits(bits::SMGNoUBBits) = bits

function smg_no_ub_all_sides(bits::SMGNoUBBits)
    return (bits.t, bits.x, bits.y, bits.tp, bits.xp, bits.yp)
end

function _smg_no_ub_pow(factor, bit)
    bit == 0 && return GrassmannExpr(1)
    return factor
end

function _smg_no_ub_leg(name)
    return grassmann(Symbol("smg_no_ub_" * name))
end

function _smg_no_ub_side_legs(prefix)
    u_forward = _smg_no_ub_leg(prefix * "_u_forward")
    u_backward = _smg_no_ub_leg(prefix * "_u_backward")
    d_forward = _smg_no_ub_leg(prefix * "_d_forward")
    d_backward = _smg_no_ub_leg(prefix * "_d_backward")
    return (; u_forward, u_backward, d_forward, d_backward)
end

function smg_no_ub_variables()
    @variables UI p_t p_x p_y

    fields = grassmann(:smg_no_ub_field, 1:4)
    physical = (;
        ubar=fields[1],
        u=fields[2],
        dbar=fields[3],
        d=fields[4],
    )
    sides = (;
        t=_smg_no_ub_side_legs("t"),
        x=_smg_no_ub_side_legs("x"),
        y=_smg_no_ub_side_legs("y"),
        tp=_smg_no_ub_side_legs("tp"),
        xp=_smg_no_ub_side_legs("xp"),
        yp=_smg_no_ub_side_legs("yp"),
    )
    parameters = (; UI, p_t, p_x, p_y)
    measure_order = GrassmannGenerator[
        physical.ubar, physical.u, physical.dbar, physical.d,
    ]
    leg_groups = (4, 4, 4, 4, 4, 4)
    return (; physical, sides, parameters, measure_order, leg_groups)
end

function smg_no_ub_onsite_fourfermi(model)
    f = model.physical
    return f.ubar * f.u * f.dbar * f.d
end

function _smg_no_ub_side_factor(
    side::SMGNoUBSideBits,
    phase,
    model;
    incoming::Bool=false,
)
    f = model.physical

    u_forward = incoming ? (-phase) * f.u : (-phase) * f.ubar
    u_backward = incoming ? phase * f.ubar : phase * f.u
    d_forward = incoming ? (-phase) * f.d : (-phase) * f.dbar
    d_backward = incoming ? phase * f.dbar : phase * f.d

    result = GrassmannExpr(1)
    result = result * _smg_no_ub_pow(u_forward, side.u_forward)
    result = result * _smg_no_ub_pow(u_backward, side.u_backward)
    result = result * _smg_no_ub_pow(d_forward, side.d_forward)
    result = result * _smg_no_ub_pow(d_backward, side.d_backward)
    return result
end

function smg_no_ub_local_integral(bits, model; simplify::Bool=true)
    local_bits = smg_no_ub_bits(bits)
    p = model.parameters
    integrand = exp(p.UI * smg_no_ub_onsite_fourfermi(model))

    for (side, phase) in zip((local_bits.t, local_bits.x, local_bits.y),
                             (p.p_t, p.p_x, p.p_y))
        integrand = integrand * _smg_no_ub_side_factor(side, phase, model)
    end
    for (side, phase) in zip((local_bits.tp, local_bits.xp, local_bits.yp),
                             (p.p_t, p.p_x, p.p_y))
        integrand = integrand *
            _smg_no_ub_side_factor(side, phase, model; incoming=true)
    end

    coefficient = scalar_part(integrate(integrand, model.measure_order))
    if simplify
        return scalar_part(simplify_coefficients(GrassmannExpr(coefficient)))
    end
    return coefficient
end

function _smg_no_ub_hopping_monomial(side::SMGNoUBSideBits, legs)
    result = GrassmannExpr(1)
    result = result *
        _smg_no_ub_pow(dual(legs.u_forward) * legs.u_forward, side.u_forward)
    result = result *
        _smg_no_ub_pow(dual(legs.u_backward) * legs.u_backward, side.u_backward)
    result = result *
        _smg_no_ub_pow(dual(legs.d_forward) * legs.d_forward, side.d_forward)
    result = result *
        _smg_no_ub_pow(dual(legs.d_backward) * legs.d_backward, side.d_backward)
    return result
end

function smg_no_ub_grassmann_monomial(bits, model)
    local_bits = smg_no_ub_bits(bits)
    result = GrassmannExpr(1)
    for (side, legs) in zip(smg_no_ub_all_sides(local_bits), values(model.sides))
        result = result * _smg_no_ub_hopping_monomial(side, legs)
    end
    return result
end

function smg_no_ub_tensor_entry(bits, model; simplify::Bool=true)
    local_bits = smg_no_ub_bits(bits)
    entry = smg_no_ub_local_integral(local_bits, model; simplify) *
        smg_no_ub_grassmann_monomial(local_bits, model)
    return simplify ? simplify_coefficients(entry) : entry
end

function derive_3d_staggered_fermion_smg_no_ub_tensor()
    model = smg_no_ub_variables()
    entry = bits -> smg_no_ub_tensor_entry(bits, model)
    coefficient = bits -> smg_no_ub_local_integral(bits, model)
    return (; model..., entry, coefficient)
end

function smg_no_ub_staggered_phase_substitutions(
    model;
    t::Integer,
    x::Integer,
    y::Integer=0,
)
    p = model.parameters
    return Dict(
        p.p_t => 1,
        p.p_x => (-1)^t,
        p.p_y => (-1)^(t + x),
    )
end

function validate_3d_staggered_fermion_smg_no_ub_tensor()
    @testset "3D staggered fermion SMG no-UB tensor" begin
        model = derive_3d_staggered_fermion_smg_no_ub_tensor()
        p = model.parameters
        @test coefficient_array_shape(model.leg_groups) == (16, 16, 16, 16, 16, 16)

        empty_bits = smg_no_ub_bits(ntuple(_ -> 0, 24))
        @test symbolically_equal(
            GrassmannExpr(model.coefficient(empty_bits)),
            GrassmannExpr(p.UI),
        )
        @test is_even(model.entry(empty_bits))

        single_u_hop = smg_no_ub_bits((1, 0, 0, 0, ntuple(_ -> 0, 20)...))
        @test model.coefficient(single_u_hop) == 0

        saturated_hops = smg_no_ub_bits((1, 1, 1, 1, ntuple(_ -> 0, 20)...))
        @test symbolically_equal(
            GrassmannExpr(model.coefficient(saturated_hops)),
            GrassmannExpr(p.p_t^4),
        )
    end
    return derive_3d_staggered_fermion_smg_no_ub_tensor()
end

if abspath(PROGRAM_FILE) == @__FILE__
    model = validate_3d_staggered_fermion_smg_no_ub_tensor()
    println("3d_staggered_fermion_smg_no_ub_shape_", coefficient_array_shape(model.leg_groups))
end
