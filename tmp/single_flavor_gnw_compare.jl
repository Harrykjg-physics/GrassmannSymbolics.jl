using Pkg
const PROJECT_ROOT = dirname(@__DIR__)
Pkg.activate(PROJECT_ROOT)
using GrassmannSymbolics
using Symbolics

include(joinpath(PROJECT_ROOT, "examples", "Single_Flavor_Gross_Neveu_Wilson.jl"))
include(joinpath(PROJECT_ROOT, "external_tensor.jl"))

function numeric_value(x)
    y = try Symbolics.simplify(x; expand=true) catch; x end
    try
        return ComplexF64(y)
    catch
        return Core.eval(Main, Meta.parse(string(y)))
    end
end

function as_complex_array(array)
    out = Array{ComplexF64}(undef, size(array))
    for I in CartesianIndices(array)
        out[I] = ComplexF64(numeric_value(array[I]))
    end
    return out
end

function arr(tensor, groups, rules)
    as_complex_array(coefficient_array(tensor, groups; substitutions=rules, element_type=Any, simplify=true))
end

function report(label, lhs, rhs)
    diff = abs.(lhs .- rhs)
    maxv = maximum(diff)
    idx = argmax(diff)
    println(label)
    println("  max_abs_diff = ", maxv)
    println("  max_index = ", Tuple(idx))
    println("  ours = ", lhs[idx])
    println("  grassmanntn = ", rhs[idx])
    println("  nonzero_diff_count_1e-10 = ", count(>(1e-10), diff))
end

model = derive_single_flavor_gnw_tensor()
D,g2,z,s = model.parameters
mu = 0.31
m = 0.37
g = 0.44
rules = Dict(D=>m+2, g2=>g, z=>exp(mu/2), s=>1/sqrt(2))
ext = Gross_Neveu_model_external_field(mu, m, g, 0.0)
legs = model.legs
println("legs = ", legs)
variants = (
    model_pair = [collect(legs[1:2]), collect(legs[3:4]), collect(legs[5:6]), collect(legs[7:8])],
    bar_forward = [collect(legs[1:2]), collect(legs[3:4]), [legs[6], legs[5]], [legs[8], legs[7]]],
)
for name in keys(variants)
    report(string(name), arr(model.tensor, getproperty(variants, name), rules), ext)
end