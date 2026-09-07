using GrassmannSymbolics
using Symbolics
using Test

# Local tensor entries for Eq. (1) of
# "Phase diagram of a lattice fermion model with symmetric mass generation".
#
# The action contains two massless staggered flavors, u and d, on a 3D cubic
# lattice.  For each incident lattice side, the six local occupation bits are
#
#   (u_forward, u_backward, u_dimer, d_forward, d_backward, d_dimer).
#
# The directions are ordered as (t, x, y, tp, xp, yp), where primed directions
# are the bonds entering the site from the negative direction.  The dimer bits
# are bosonic occupation indices for the UB four-fermi bond term.  The local
# endpoint weight is sqrt_UB * (bar(f) f), so contracting the same dimer bit
# across the two bond endpoints reconstructs UB * (bar(f) f)_i (bar(f) f)_j.

struct SMGSideBits
    u_forward::Int
    u_backward::Int
    u_dimer::Int
    d_forward::Int
    d_backward::Int
    d_dimer::Int
end

struct SMGBits
    t::SMGSideBits
    x::SMGSideBits
    y::SMGSideBits
    tp::SMGSideBits
    xp::SMGSideBits
    yp::SMGSideBits
end

function _smg_check_bit(bit)
    bit in (0, 1) || throw(ArgumentError("SMG occupation bits must be 0 or 1"))
    return Int(bit)
end

function smg_side_bits(values)
    length(values) == 6 ||
        throw(ArgumentError("each SMG side needs 6 occupation bits"))
    checked = Tuple(_smg_check_bit(bit) for bit in values)
    return SMGSideBits(checked...)
end

function smg_bits(values)
    length(values) == 36 ||
        throw(ArgumentError("SMG local tensor entries need 36 occupation bits"))
    return SMGBits(
        smg_side_bits(values[1:6]),
        smg_side_bits(values[7:12]),
        smg_side_bits(values[13:18]),
        smg_side_bits(values[19:24]),
        smg_side_bits(values[25:30]),
        smg_side_bits(values[31:36]),
    )
end

smg_bits(bits::SMGBits) = bits

smg_all_sides(bits::SMGBits) = (bits.t, bits.x, bits.y, bits.tp, bits.xp, bits.yp)

function _smg_pow(factor, bit)
    bit == 0 && return GrassmannExpr(1)
    return factor
end

function _smg_leg(name)
    return grassmann(Symbol("smg_" * name))
end

function _smg_side_legs(prefix)
    u_forward = _smg_leg(prefix * "_u_forward")
    u_backward = _smg_leg(prefix * "_u_backward")
    d_forward = _smg_leg(prefix * "_d_forward")
    d_backward = _smg_leg(prefix * "_d_backward")
    return (; u_forward, u_backward, d_forward, d_backward)
end

function smg_variables()
    @variables UI sqrt_UB p_t p_x p_y

    fields = grassmann(:smg_field, 1:4)
    physical = (;
        ubar=fields[1],
        u=fields[2],
        dbar=fields[3],
        d=fields[4],
    )
    sides = (;
        t=_smg_side_legs("t"),
        x=_smg_side_legs("x"),
        y=_smg_side_legs("y"),
        tp=_smg_side_legs("tp"),
        xp=_smg_side_legs("xp"),
        yp=_smg_side_legs("yp"),
    )
    parameters = (; UI, sqrt_UB, p_t, p_x, p_y)
    measure_order = GrassmannGenerator[
        physical.ubar, physical.u, physical.dbar, physical.d,
    ]
    leg_groups = (6, 6, 6, 6, 6, 6)
    return (; physical, sides, parameters, measure_order, leg_groups)
end

smg_number_u(model) = model.physical.ubar * model.physical.u
smg_number_d(model) = model.physical.dbar * model.physical.d
smg_onsite_fourfermi(model) = smg_number_u(model) * smg_number_d(model)

function _smg_side_factor(side::SMGSideBits, phase, model; incoming::Bool=false)
    p = model.parameters
    f = model.physical

    u_forward = incoming ? (-phase) * f.u : (-phase) * f.ubar
    u_backward = incoming ? phase * f.ubar : phase * f.u
    d_forward = incoming ? (-phase) * f.d : (-phase) * f.dbar
    d_backward = incoming ? phase * f.dbar : phase * f.d

    result = GrassmannExpr(1)
    result = result * _smg_pow(u_forward, side.u_forward)
    result = result * _smg_pow(u_backward, side.u_backward)
    result = result * _smg_pow(p.sqrt_UB * smg_number_u(model), side.u_dimer)
    result = result * _smg_pow(d_forward, side.d_forward)
    result = result * _smg_pow(d_backward, side.d_backward)
    result = result * _smg_pow(p.sqrt_UB * smg_number_d(model), side.d_dimer)
    return result
end

function smg_local_integral(bits, model; simplify::Bool=true)
    local_bits = smg_bits(bits)
    p = model.parameters
    integrand = exp(p.UI * smg_onsite_fourfermi(model))

    for (side, phase) in zip((local_bits.t, local_bits.x, local_bits.y),
                             (p.p_t, p.p_x, p.p_y))
        integrand = integrand * _smg_side_factor(side, phase, model)
    end
    for (side, phase) in zip((local_bits.tp, local_bits.xp, local_bits.yp),
                             (p.p_t, p.p_x, p.p_y))
        integrand = integrand * _smg_side_factor(side, phase, model; incoming=true)
    end

    coefficient = scalar_part(integrate(integrand, model.measure_order))
    if simplify
        return scalar_part(simplify_coefficients(GrassmannExpr(coefficient)))
    end
    return coefficient
end

function _smg_hopping_monomial(side::SMGSideBits, legs)
    result = GrassmannExpr(1)
    result = result * _smg_pow(dual(legs.u_forward) * legs.u_forward, side.u_forward)
    result = result * _smg_pow(dual(legs.u_backward) * legs.u_backward, side.u_backward)
    result = result * _smg_pow(dual(legs.d_forward) * legs.d_forward, side.d_forward)
    result = result * _smg_pow(dual(legs.d_backward) * legs.d_backward, side.d_backward)
    return result
end

function smg_grassmann_monomial(bits, model)
    local_bits = smg_bits(bits)
    result = GrassmannExpr(1)
    for (side, legs) in zip(smg_all_sides(local_bits), values(model.sides))
        result = result * _smg_hopping_monomial(side, legs)
    end
    return result
end

function smg_tensor_entry(bits, model; simplify::Bool=true)
    local_bits = smg_bits(bits)
    entry = smg_local_integral(local_bits, model; simplify) *
        smg_grassmann_monomial(local_bits, model)
    return simplify ? simplify_coefficients(entry) : entry
end

function derive_3d_staggered_fermion_smg_tensor()
    model = smg_variables()
    entry = bits -> smg_tensor_entry(bits, model)
    coefficient = bits -> smg_local_integral(bits, model)
    return (; model..., entry, coefficient)
end

function smg_staggered_phase_substitutions(model; t::Integer, x::Integer, y::Integer=0)
    p = model.parameters
    return Dict(
        p.p_t => 1,
        p.p_x => (-1)^t,
        p.p_y => (-1)^(t + x),
    )
end

function validate_3d_staggered_fermion_smg_tensor()
    @testset "3D staggered fermion SMG tensor" begin
        model = derive_3d_staggered_fermion_smg_tensor()
        p = model.parameters
        @test coefficient_array_shape(model.leg_groups) == (64, 64, 64, 64, 64, 64)

        empty_bits = smg_bits(ntuple(_ -> 0, 36))
        @test symbolically_equal(
            GrassmannExpr(model.coefficient(empty_bits)),
            GrassmannExpr(p.UI),
        )
        @test is_even(model.entry(empty_bits))

        single_u_dimer = smg_bits((0, 0, 1, 0, 0, 0, ntuple(_ -> 0, 30)...))
        @test model.coefficient(single_u_dimer) == 0

        complementary_dimers = smg_bits((0, 0, 1, 0, 0, 1, ntuple(_ -> 0, 30)...))
        @test symbolically_equal(
            GrassmannExpr(model.coefficient(complementary_dimers)),
            GrassmannExpr(p.sqrt_UB^2),
        )

        repeated_u_dimer = smg_bits((0, 0, 1, 0, 0, 0,
                                     0, 0, 0, 0, 0, 0,
                                     0, 0, 0, 0, 0, 0,
                                     0, 0, 1, 0, 0, 0,
                                     0, 0, 0, 0, 0, 0,
                                     0, 0, 0, 0, 0, 0))
        @test model.coefficient(repeated_u_dimer) == 0
    end
    return derive_3d_staggered_fermion_smg_tensor()
end

if abspath(PROGRAM_FILE) == @__FILE__
    model = validate_3d_staggered_fermion_smg_tensor()
    println("3d_staggered_fermion_smg_shape_", coefficient_array_shape(model.leg_groups))
end
