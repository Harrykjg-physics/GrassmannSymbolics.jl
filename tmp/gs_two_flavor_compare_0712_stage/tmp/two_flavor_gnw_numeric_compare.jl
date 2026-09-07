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
    (; label="mu_zero", 渭1=0.0, 渭2=0.0, m1=0.37, m2=-0.21, gs2=0.44, gp2=0.13, H=0.17),
    (; label="mu_nonzero", 渭1=0.31, 渭2=-0.23, m1=0.37, m2=-0.21, gs2=0.44, gp2=0.13, H=0.17),
]

function tfgnw_groups(model)
    return [
        collect(model.legs[1:4]),
        collect(model.legs[5:8]),
        collect(model.legs[9:12]),
        collect(model.legs[13:16]),
    ]
end

function tfgnw_rules(model, params)
    p = model.params
    return Dict(
        p.M1 => params.m1,
        p.M2 => params.m2,
        p.gs2 => params.gs2,
        p.gp2 => params.gp2,
        p.H => params.H,
        p.q => 1 / sqrt(2),
        p.z1 => exp(params.渭1 / 2),
        p.z2 => exp(params.渭2 / 2),
    )
end

function as_complex_array(array)
    out = Array{ComplexF64}(undef, size(array))
    for index in CartesianIndices(array)
        out[index] = ComplexF64(array[index])
    end
    return out
end

function derived_arrays(model, params)
    groups = tfgnw_groups(model)
    rules = tfgnw_rules(model, params)
    return (;
        pure = coefficient_array(model.pure, groups; substitutions=rules, element_type=ComplexF64),
        chiral = coefficient_array(model.chiral_condensate, groups; substitutions=rules, element_type=ComplexF64),
        singlet = coefficient_array(model.pseudoscalar_singlet, groups; substitutions=rules, element_type=ComplexF64),
        triplet = coefficient_array(model.pseudoscalar_triplet, groups; substitutions=rules, element_type=ComplexF64),
    )
end

function maxdiff_report(label, lhs, rhs)
    diff = abs.(lhs .- rhs)
    max_value = maximum(diff)
    max_index = argmax(diff)
    println(label)
    println("  max_abs_diff = ", max_value)
    println("  max_index = ", Tuple(max_index))
    println("  ours = ", lhs[max_index])
    println("  grassmanntn = ", rhs[max_index])
    println("  nonzero_diff_count_1e-10 = ", count(>(1e-10), diff))
    return max_value
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
    ours = derived_arrays(model, params)
    external_pure = Gross_Neveu_Wilson_model_Nf2(
        params.渭1, params.渭2, params.m1, params.m2, params.gs2, params.gp2, params.H,
    )
    maxdiff_report("pure", ours.pure, external_pure)

    external_chiral, external_singlet, external_triplet =
        Gross_Neveu_Wilson_model_Nf2_impurity(
            params.m1, params.m2, params.gs2, params.gp2, params.H,
        )
    maxdiff_report("chiral_impurity", ours.chiral, external_chiral)
    maxdiff_report("pseudoscalar_singlet_impurity", ours.singlet, external_singlet)
    maxdiff_report("pseudoscalar_triplet_impurity", ours.triplet, external_triplet)
end




