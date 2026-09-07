using GrassmannSymbolics
using Symbolics
using Test

# Three-flavor Gross-Neveu-Wilson model from Ref/2d_gn.pdf, section 4.
#
# Nf=3 has 24 external Grassmann legs.  A fully expanded sparse GrassmannExpr
# would require scanning a 2^24-dimensional occupation space, so the practical
# tensor representation here is a symbolic coefficient oracle: for a given
# 24-bit external occupation pattern it selects the corresponding Taylor factors
# in Eq. (4.4) and performs the local Berezin integral exactly.  This is the
# same representation used by the long coefficient tables in the reference.
#
# The two impurity insertions in Eqs. (4.68)--(4.70) are implemented as
# observable=:pseudoscalar_singlet_squared and
# observable=:pseudoscalar_lambda8_squared.

struct ThreeFlavorGNWModel
    vars::Any
    params::Any
    onsite_weight::GrassmannExpr
    leg_exponents::Vector{Pair{GrassmannGenerator,Any}}
    insertions::Any
    legs::Vector{GrassmannGenerator}
    measure_order::Vector{GrassmannGenerator}
end

function _thfgnw_generators()
    psi = [grassmann(Symbol(:psi, f, s)) for f in 1:3, s in 1:2]
    psibar = dual.(psi)

    alpha = grassmann(:alpha, 1:3)
    beta = grassmann(:beta, 1:3)
    eta = grassmann(:eta, 1:3)
    zeta = grassmann(:zeta, 1:3)
    alphabar = dual.(alpha)
    betabar = dual.(beta)
    etabar = dual.(eta)
    zetabar = dual.(zeta)

    measure_order = GrassmannGenerator[]
    for f in 1:3
        append!(measure_order, (psi[f, 1], psibar[f, 1], psi[f, 2], psibar[f, 2]))
    end

    # Flattened convention parallel to Eq. (4.1): flavor blocks first, then the
    # dual legs in the order used for contractions.
    legs = GrassmannGenerator[]
    for f in 1:3
        append!(legs, (alpha[f], beta[f], eta[f], zeta[f]))
    end
    for f in 1:3
        append!(legs, (zetabar[f], etabar[f], betabar[f], alphabar[f]))
    end

    return (; psi, psibar, alpha, beta, eta, zeta, alphabar, betabar,
        etabar, zetabar, measure_order, legs)
end

function _thfgnw_flavor_fields(vars, f, q)
    psi = vars.psi
    psibar = vars.psibar

    phi1 = q * (psi[f, 1] + psi[f, 2])
    phi2 = q * (psi[f, 1] - psi[f, 2])
    phibar1 = q * (psibar[f, 1] + psibar[f, 2])
    phibar2 = q * (psibar[f, 1] - psibar[f, 2])

    chi1 = q * (psi[f, 1] - im * psi[f, 2])
    chi2 = q * (psi[f, 1] + im * psi[f, 2])
    chibar1 = q * (psibar[f, 1] + im * psibar[f, 2])
    chibar2 = q * (psibar[f, 1] - im * psibar[f, 2])

    return (; phi1, phi2, phibar1, phibar2, chi1, chi2, chibar1, chibar2)
end

_thfgnw_density(vars, f, s) = vars.psibar[f, s] * vars.psi[f, s]

function _thfgnw_onsite_action(vars, params)
    masses = (params.M1, params.M2, params.M3)
    axial_charges = (1, 1, -2)
    H, gs2, gp2 = params.H, params.gs2, params.gp2

    densities = Dict((f, s) => _thfgnw_density(vars, f, s) for f in 1:3 for s in 1:2)
    mass_action = GrassmannExpr(0)
    for f in 1:3
        charge = axial_charges[f]
        mass_action = mass_action +
            (masses[f] + 2 + im * charge * H) * densities[(f, 1)] +
            (masses[f] + 2 - im * charge * H) * densities[(f, 2)]
    end

    equal_chirality = GrassmannExpr(0)
    opposite_chirality = GrassmannExpr(0)
    for f in 1:3
        opposite_chirality = opposite_chirality + densities[(f, 1)] * densities[(f, 2)]
    end
    for f in 1:3, g in (f + 1):3
        equal_chirality = equal_chirality +
            densities[(f, 1)] * densities[(g, 1)] +
            densities[(f, 2)] * densities[(g, 2)]
        opposite_chirality = opposite_chirality +
            densities[(f, 1)] * densities[(g, 2)] +
            densities[(f, 2)] * densities[(g, 1)]
    end

    interaction_exponent = ((1 // 3) * (gs2 - gp2)) * equal_chirality +
        ((1 // 3) * (gs2 + gp2)) * opposite_chirality
    return mass_action - interaction_exponent
end

function _thfgnw_onsite_weight(vars, params)
    masses = (params.M1, params.M2, params.M3)
    axial_charges = (1, 1, -2)
    H, gs2, gp2 = params.H, params.gs2, params.gp2
    equal_coeff = (1 // 3) * (gs2 - gp2)
    opposite_coeff = (1 // 3) * (gs2 + gp2)

    densities = Dict((f, s) => _thfgnw_density(vars, f, s) for f in 1:3 for s in 1:2)
    weight = GrassmannExpr(1)

    # Exact factorized form of exp(-S_onsite).  The density monomials are even
    # and nilpotent, so exp(c*n)=1+c*n and exp(c*n_i*n_j)=1+c*n_i*n_j.
    for f in 1:3
        charge = axial_charges[f]
        a1 = masses[f] + 2 + im * charge * H
        a2 = masses[f] + 2 - im * charge * H
        weight = weight * (GrassmannExpr(1) - a1 * densities[(f, 1)])
        weight = weight * (GrassmannExpr(1) - a2 * densities[(f, 2)])
    end

    for f in 1:3, g in (f + 1):3
        weight = weight * (GrassmannExpr(1) + equal_coeff * densities[(f, 1)] * densities[(g, 1)])
        weight = weight * (GrassmannExpr(1) + equal_coeff * densities[(f, 2)] * densities[(g, 2)])
        weight = weight * (GrassmannExpr(1) + opposite_coeff * densities[(f, 1)] * densities[(g, 2)])
        weight = weight * (GrassmannExpr(1) + opposite_coeff * densities[(f, 2)] * densities[(g, 1)])
    end
    for f in 1:3
        weight = weight * (GrassmannExpr(1) + opposite_coeff * densities[(f, 1)] * densities[(f, 2)])
    end

    return weight
end

function _thfgnw_leg_exponents(vars, params)
    q = params.q
    z = (params.z1, params.z2, params.z3)
    pairs = Pair{GrassmannGenerator,Any}[]

    # Product order follows the local factors in Eq. (4.4).  The final tensor
    # leg order is stored separately in model.legs.
    for f in 1:3
        fields = _thfgnw_flavor_fields(vars, f, q)
        push!(pairs, vars.alpha[f] => fields.phibar2 * vars.alpha[f])
        push!(pairs, vars.alphabar[f] => -fields.phi2 * vars.alphabar[f])
        push!(pairs, vars.betabar[f] => fields.phibar1 * vars.betabar[f])
        push!(pairs, vars.beta[f] => fields.phi1 * vars.beta[f])

        push!(pairs, vars.eta[f] => z[f] * fields.chibar2 * vars.eta[f])
        push!(pairs, vars.etabar[f] => -z[f] * fields.chi2 * vars.etabar[f])
        push!(pairs, vars.zetabar[f] => (1 / z[f]) * fields.chibar1 * vars.zetabar[f])
        push!(pairs, vars.zeta[f] => (1 / z[f]) * fields.chi1 * vars.zeta[f])
    end

    return pairs
end

function _thfgnw_insertions(vars)
    p = [_thfgnw_density(vars, f, 1) - _thfgnw_density(vars, f, 2) for f in 1:3]
    pseudoscalar_singlet_squared = -(p[1] + p[2] + p[3])^2
    pseudoscalar_lambda8_squared = -((1 // 3) * (p[1] + p[2] - 2 * p[3])^2)
    return (; pseudoscalar_singlet_squared, pseudoscalar_lambda8_squared)
end

function three_flavor_gnw_model()
    @variables M1 M2 M3 gs2 gp2 H q z1 z2 z3
    params = (; M1, M2, M3, gs2, gp2, H, q, z1, z2, z3)
    vars = _thfgnw_generators()
    onsite_weight = _thfgnw_onsite_weight(vars, params)
    leg_exponents = _thfgnw_leg_exponents(vars, params)
    insertions = _thfgnw_insertions(vars)
    return ThreeFlavorGNWModel(
        vars,
        params,
        onsite_weight,
        leg_exponents,
        insertions,
        vars.legs,
        vars.measure_order,
    )
end

function _thfgnw_observable_factor(model::ThreeFlavorGNWModel, observable::Symbol)
    observable === :pure && return GrassmannExpr(1)
    observable === :pseudoscalar_singlet_squared &&
        return model.insertions.pseudoscalar_singlet_squared
    observable === :pseudoscalar_lambda8_squared &&
        return model.insertions.pseudoscalar_lambda8_squared
    throw(ArgumentError("unknown three-flavor GNW observable: $(observable)"))
end

function three_flavor_gnw_coefficient(
    model::ThreeFlavorGNWModel,
    bits;
    observable::Symbol=:pure,
    simplify::Bool=true,
)
    length(bits) == 24 || throw(ArgumentError("three-flavor GNW expects 24 leg bits"))
    occupied_legs = GrassmannGenerator[leg for (leg, bit) in zip(model.legs, bits) if bit == 1]
    occupied = Set(occupied_legs)

    integrand = _thfgnw_observable_factor(model, observable) * model.onsite_weight
    for (leg, exponent) in model.leg_exponents
        leg in occupied && (integrand = integrand * exponent)
    end

    projected = integrate(integrand, model.measure_order)
    coefficient = get_coeff(projected, occupied_legs)
    return simplify ? simplify_coefficients(GrassmannExpr(coefficient)) : GrassmannExpr(coefficient)
end

function three_flavor_gnw_sample_tensor(
    model::ThreeFlavorGNWModel,
    masks;
    observable::Symbol=:pure,
    simplify::Bool=true,
)
    result = Dict{NTuple{24,Int},GrassmannExpr}()
    for mask in masks
        bits = ntuple(k -> (mask >> (k - 1)) & 1, 24)
        result[bits] = three_flavor_gnw_coefficient(
            model,
            bits;
            observable=observable,
            simplify=simplify,
        )
    end
    return result
end

function validate_three_flavor_gnw_coefficients(; sample_masks=(0, 1, 3, 15, 255))
    model = three_flavor_gnw_model()
    @test length(model.legs) == 24
    @test length(model.measure_order) == 12

    pure_samples = three_flavor_gnw_sample_tensor(model, sample_masks; observable=:pure, simplify=false)
    singlet_samples = three_flavor_gnw_sample_tensor(model, sample_masks; observable=:pseudoscalar_singlet_squared, simplify=false)
    lambda8_samples = three_flavor_gnw_sample_tensor(model, sample_masks; observable=:pseudoscalar_lambda8_squared, simplify=false)

    @test haskey(pure_samples, ntuple(_ -> 0, 24))
    @test all(value -> is_even(value), values(pure_samples))
    @test all(value -> is_even(value), values(singlet_samples))
    @test all(value -> is_even(value), values(lambda8_samples))

    # The fully empty external leg component is a useful normalization check:
    # with no auxiliary legs selected, only saturated onsite monomials survive.
    empty_bits = ntuple(_ -> 0, 24)
    @test !iszero(pure_samples[empty_bits])

    return (; model, pure_samples, singlet_samples, lambda8_samples)
end

if abspath(PROGRAM_FILE) == @__FILE__
    validation = validate_three_flavor_gnw_coefficients()
    println("three_flavor_GNW_coefficient_validation_passed")
    println("sample_count_", length(validation.pure_samples))
    println("empty_pure_terms_", length(validation.pure_samples[ntuple(_ -> 0, 24)].terms))
end
