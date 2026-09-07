using Pkg
const PROJECT_ROOT = dirname(@__DIR__)
Pkg.activate(PROJECT_ROOT)
Pkg.instantiate()

using GrassmannSymbolics
using Symbolics

include(joinpath(PROJECT_ROOT, "examples", "two_flavor_Gross_Neveu_Wilson.jl"))
external_tensor_path = get(ENV, "GRASSMANNTN_TENSOR_JL", joinpath(PROJECT_ROOT, "external_tensor.jl"))
isfile(external_tensor_path) || error("external tensor.jl not found: ", external_tensor_path)
include(external_tensor_path)

const PARAM_SETS = [
    (; label="mu_zero", mu1=0.0, mu2=0.0, m1=0.37, m2=-0.21, gs2=0.44, gp2=0.13, H=0.17),
    (; label="mu_nonzero", mu1=0.31, mu2=-0.23, m1=0.37, m2=-0.21, gs2=0.44, gp2=0.13, H=0.17),
]

const GRASSMANNTN_F4_BITS = [
    (0,0,0,0), (1,1,0,0), (1,0,1,0), (1,0,0,1),
    (0,1,1,0), (0,1,0,1), (0,0,1,1), (1,1,1,1),
    (1,0,0,0), (0,1,0,0), (0,0,1,0), (1,1,1,0),
    (0,0,0,1), (1,1,0,1), (1,0,1,1), (0,1,1,1),
]

function tfgnw_rules(model, params)
    p = model.params
    return Dict(
        p.M1 => params.m1,
        p.M2 => params.m2,
        p.gs2 => params.gs2,
        p.gp2 => params.gp2,
        p.H => params.H,
        p.q => 1 / sqrt(2),
        p.z1 => exp(params.mu1 / 2),
        p.z2 => exp(params.mu2 / 2),
    )
end

function numeric_value(x)
    y = try Symbolics.simplify(x; expand=true) catch; x end
    try
        return ComplexF64(y)
    catch
        return Core.eval(Main, Meta.parse(string(y)))
    end
end

function external_index_occupied_legs(vars, index_tuple)
    b1 = GRASSMANNTN_F4_BITS[index_tuple[1]] # i1,j1,k1,l1: x direction, flavor 1 then flavor 2
    b2 = GRASSMANNTN_F4_BITS[index_tuple[2]] # i2,j2,k2,l2: y direction, flavor 1 then flavor 2
    b3 = GRASSMANNTN_F4_BITS[index_tuple[3]] # i1p,j1p,k1p,l1p: x-bar direction
    b4 = GRASSMANNTN_F4_BITS[index_tuple[4]] # i2p,j2p,k2p,l2p: y-bar direction
    i1,j1,k1,l1 = b1
    i2,j2,k2,l2 = b2
    i1p,j1p,k1p,l1p = b3
    i2p,j2p,k2p,l2p = b4

    ordered = GrassmannGenerator[
        vars.alpha[1], vars.beta[1], vars.alpha[2], vars.beta[2],
        vars.eta[1], vars.zeta[1], vars.eta[2], vars.zeta[2],
        vars.betabar[1], vars.alphabar[1], vars.betabar[2], vars.alphabar[2],
        vars.zetabar[1], vars.etabar[1], vars.zetabar[2], vars.etabar[2],
    ]
    occupations = (i1,j1,k1,l1, i2,j2,k2,l2, j1p,i1p,l1p,k1p, j2p,i2p,l2p,k2p)
    return GrassmannGenerator[g for (g, bit) in zip(ordered, occupations) if bit == 1]
end

function tensor_array_external_basis(tensor, model, rules)
    prepared = substitute_coefficients(tensor, rules)
    out = Array{ComplexF64}(undef, 16, 16, 16, 16)
    for I in CartesianIndices(out)
        coeff = get_coeff(prepared, external_index_occupied_legs(model.vars, Tuple(I)))
        out[I] = ComplexF64(numeric_value(coeff))
    end
    return out
end

function print_report(label, lhs, rhs)
    diff = abs.(lhs .- rhs)
    max_value = maximum(diff)
    max_index = argmax(diff)
    println(label)
    println("  max_abs_diff = ", max_value)
    println("  max_index = ", Tuple(max_index))
    println("  ours = ", lhs[max_index])
    println("  grassmanntn = ", rhs[max_index])
    println("  nonzero_diff_count_1e-10 = ", count(>(1e-10), diff))
end

println("derive_two_flavor_gnw_tensors started")
model = derive_two_flavor_gnw_tensors()
println("derive_two_flavor_gnw_tensors done")
println("terms pure=", length(model.pure.terms),
    " chiral=", length(model.chiral_condensate.terms),
    " singlet=", length(model.pseudoscalar_singlet.terms),
    " triplet=", length(model.pseudoscalar_triplet.terms))

for params in PARAM_SETS
    println()
    println("PARAM_SET ", params.label)
    println(params)
    rules = tfgnw_rules(model, params)
    external_pure = Gross_Neveu_Wilson_model_Nf2(
        params.mu1, params.mu2, params.m1, params.m2, params.gs2, params.gp2, params.H,
    )
    print_report("pure_direct_external_basis", tensor_array_external_basis(model.pure, model, rules), external_pure)
end

params = PARAM_SETS[1]
println()
println("IMPURITY_AT_MU_ZERO")
println(params)
rules = tfgnw_rules(model, params)
external_chiral, external_singlet, external_triplet = Gross_Neveu_Wilson_model_Nf2_impurity(
    params.m1, params.m2, params.gs2, params.gp2, params.H,
)
print_report("chiral_direct_external_basis", tensor_array_external_basis(model.chiral_condensate, model, rules), external_chiral)
print_report("pseudoscalar_singlet_direct_external_basis", tensor_array_external_basis(model.pseudoscalar_singlet, model, rules), external_singlet)
print_report("pseudoscalar_triplet_direct_external_basis", tensor_array_external_basis(model.pseudoscalar_triplet, model, rules), external_triplet)