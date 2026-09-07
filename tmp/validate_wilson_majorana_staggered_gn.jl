using Pkg
Pkg.activate(".")
Pkg.instantiate()

include(joinpath(pwd(), "examples", "Wilson_Majorana_Fermion.jl"))
include(joinpath(pwd(), "examples", "Two_Flavor_Staggered_Gross_Neveu.jl"))
include(joinpath(pwd(), "examples", "coefficient_tensors.jl"))

wm = validate_wilson_majorana_tensor()
tf = validate_two_flavor_staggered_gn_tensor()

wm_spec = wilson_majorana_coefficient_tensor(; materialize=false)
tf_spec = two_flavor_staggered_gn_coefficient_tensor(; materialize=false)

@assert wm_spec.shape == (16, 16, 16, 16)
@assert tf_spec.shape == (16, 16, 16, 16)
@assert length(wm.tensor.terms) > 0
@assert length(tf.tensor.terms) > 0

println("WILSON_MAJORANA_TERMS=", length(wm.tensor.terms))
println("TWO_FLAVOR_STAGGERED_GN_TERMS=", length(tf.tensor.terms))
println("WM_SPEC_SHAPE=", wm_spec.shape)
println("TFSGN_SPEC_SHAPE=", tf_spec.shape)

Pkg.activate("docs")
Pkg.develop(PackageSpec(path=pwd()))
Pkg.instantiate()
include(joinpath(pwd(), "docs", "make.jl"))
println("DOCS_BUILD_DONE")
