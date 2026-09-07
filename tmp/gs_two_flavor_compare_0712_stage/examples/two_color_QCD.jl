using GrassmannSymbolics
using Symbolics
using Test

# Appendix A of JHEP 03 (2025) 027 for the N = 2 staggered fermion sector.
# The four tensor legs are X, T, Xbar, Tbar; each leg is a 4-bit color/source
# super-index.  The symbols ap and am denote a_+ and a_- respectively.

function two_color_qcd_variables()
    @variables m λ ap am
    @variables Ux11 Ux12 Ux21 Ux22 Ut11 Ut12 Ut21 Ut22
    @variables Uxd11 Uxd12 Uxd21 Uxd22 Utd11 Utd12 Utd21 Utd22

    χ1 = grassmann(:χ, index=1)
    χ2 = grassmann(:χ, index=2)
    χbar1 = dual(χ1)
    χbar2 = dual(χ2)

    ηx1 = grassmann(:ηx, index=1)
    ηx2 = grassmann(:ηx, index=2)
    ζx1 = grassmann(:ζx, index=1)
    ζx2 = grassmann(:ζx, index=2)
    ηt1 = grassmann(:ηt, index=1)
    ηt2 = grassmann(:ηt, index=2)
    ζt1 = grassmann(:ζt, index=1)
    ζt2 = grassmann(:ζt, index=2)

    ηbarx1 = dual(ηx1)
    ηbarx2 = dual(ηx2)
    ζbarx1 = dual(ζx1)
    ζbarx2 = dual(ζx2)
    ηbart1 = dual(ηt1)
    ηbart2 = dual(ηt2)
    ζbart1 = dual(ζt1)
    ζbart2 = dual(ζt2)

    original = (; χ1, χbar1, χ2, χbar2)
    auxiliaries = (;
        ηx1, ηx2, ζx1, ζx2, ηt1, ηt2, ζt1, ζt2,
        ζbarx2, ζbarx1, ηbarx2, ηbarx1,
        ζbart2, ζbart1, ηbart2, ηbart1,
    )
    gauge = (;
        Ux11, Ux12, Ux21, Ux22, Ut11, Ut12, Ut21, Ut22,
        Uxd11, Uxd12, Uxd21, Uxd22, Utd11, Utd12, Utd21, Utd22,
    )
    parameters = (; m, λ, ap, am, gauge...)
    return (; original, auxiliaries, parameters)
end

function two_color_qcd_appendix_ABCD(vars)
    a = vars.auxiliaries
    p = vars.parameters

    A = (
        -a.ηx1 - a.ηt1 +
        p.Uxd11 * a.ζbarx1 + p.Uxd12 * a.ζbarx2 +
        p.Utd11 * a.ζbart1 + p.Utd12 * a.ζbart2
    )
    B = (1 // 2) * (
        -p.Ux11 * a.ηbarx1 - p.Ux21 * a.ηbarx2 -
        p.ap * p.Ut11 * a.ηbart1 - p.ap * p.Ut21 * a.ηbart2 +
        a.ζx1 + p.am * a.ζt1
    )
    C = (
        -a.ηx2 - a.ηt2 +
        p.Uxd21 * a.ζbarx1 + p.Uxd22 * a.ζbarx2 +
        p.Utd21 * a.ζbart1 + p.Utd22 * a.ζbart2
    )
    D = (1 // 2) * (
        -p.Ux12 * a.ηbarx1 - p.Ux22 * a.ηbarx2 -
        p.ap * p.Ut12 * a.ηbart1 - p.ap * p.Ut22 * a.ηbart2 +
        a.ζx2 + p.am * a.ζt2
    )
    return (; A, B, C, D)
end

function derive_two_color_qcd_fermion_tensor(; diquark_source=false)
    vars = two_color_qcd_variables()
    o = vars.original
    p = vars.parameters
    abcd = two_color_qcd_appendix_ABCD(vars)

    onsite_action = p.m * (o.χbar1 * o.χ1 + o.χbar2 * o.χ2)
    if diquark_source
        onsite_action = onsite_action - im * p.λ * (o.χ1 * o.χ2 + o.χbar1 * o.χbar2)
    end

    bond_exponents = [
        o.χbar1 * abcd.A,
        o.χ1 * abcd.B,
        o.χbar2 * abcd.C,
        o.χ2 * abcd.D,
    ]
    tensor = local_grassmann_tensor(
        onsite_action,
        bond_exponents,
        [o.χ1, o.χbar1, o.χ2, o.χbar2],
    )

    a = vars.auxiliaries
    legs = (
        a.ηx1, a.ηx2, a.ζx1, a.ζx2,
        a.ηt1, a.ηt2, a.ζt1, a.ζt2,
        a.ζbarx2, a.ζbarx1, a.ηbarx2, a.ηbarx1,
        a.ζbart2, a.ζbart1, a.ηbart2, a.ηbart1,
    )
    return (; tensor, legs, vars, abcd, diquark_source)
end

function published_two_color_qcd_fermion_tensor(model)
    A, B, C, D = model.abcd
    m = model.vars.parameters.m
    λ = model.vars.parameters.λ
    result = A * B * C * D + m * (A * B + C * D) + m^2
    if model.diquark_source
        result = result + im * λ * (A * C + B * D) + λ^2
    end
    return simplify_coefficients(result)
end

function two_color_qcd_occupation_bits(mask)
    return ntuple(k -> (mask >> (k - 1)) & 1, 16)
end

function two_color_qcd_occupied_legs(legs, bits)
    return GrassmannGenerator[
        leg for (leg, occupation) in zip(legs, bits) if occupation == 1
    ]
end

# Equation (2.13): the bosonic plaquette weight used when beta is finite.
function two_color_qcd_gauge_plaquette_factor(Ui, Uj, Uk, Ul, β, K; N=2)
    return exp((β / N) * real(tr(Ui * adjoint(Uj) * adjoint(Uk) * Ul))) / K^2
end

# Equation (2.15): q = (q_f, q_g), so a full initial entry is F times G.
function two_color_qcd_initial_tensor_entry(fermion_entry, plaquette_entry)
    return fermion_entry * plaquette_entry
end

function validate_two_color_qcd_fermion_tensor()
    model = derive_two_color_qcd_fermion_tensor()
    reference = published_two_color_qcd_fermion_tensor(model)
    @test symbolically_equal(model.tensor, reference)
    @test is_even(model.tensor)
    @test symbolically_equal(
        GrassmannExpr(get_coeff(model.tensor, GrassmannGenerator[])),
        GrassmannExpr(model.vars.parameters.m^2),
    )

    sourced = derive_two_color_qcd_fermion_tensor(; diquark_source=true)
    sourced_reference = published_two_color_qcd_fermion_tensor(sourced)
    @test symbolically_equal(sourced.tensor, sourced_reference)
    @test is_even(sourced.tensor)

    for mask in (0, 1, 3, 15, 255, 4095, 65535)
        bits = two_color_qcd_occupation_bits(mask)
        occupied = two_color_qcd_occupied_legs(model.legs, bits)
        actual = get_coeff(model.tensor, occupied)
        expected = get_coeff(reference, occupied)
        @test symbolically_equal(GrassmannExpr(actual), GrassmannExpr(expected))
    end
    return model
end

if abspath(PROGRAM_FILE) == @__FILE__
    validate_two_color_qcd_fermion_tensor()
end
