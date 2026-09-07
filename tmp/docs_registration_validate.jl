using Pkg

const PROJECT_ROOT = pwd()

println("JULIA_VERSION=", VERSION)
println("PROJECT=", PROJECT_ROOT)

println("== package instantiate ==")
Pkg.instantiate()

println("== package test ==")
Pkg.test()

println("== docs instantiate ==")
Pkg.activate("docs")
Pkg.develop(PackageSpec(path=PROJECT_ROOT))
Pkg.instantiate()

println("== docs build ==")
include(joinpath(PROJECT_ROOT, "docs", "make.jl"))

println("VALIDATION_DONE")
