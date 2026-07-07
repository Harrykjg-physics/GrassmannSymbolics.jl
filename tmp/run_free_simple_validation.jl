import Pkg

Pkg.activate(@__DIR__)

include(joinpath(@__DIR__, string(:test), string(:runtests, Char(46), :jl)))
include(joinpath(
    @__DIR__,
    string(:examples),
    string(:Free_Wilson_and_Staggered, Char(46), :jl),
))
include(joinpath(
    @__DIR__,
    string(:examples),
    string(:Simple_quardratic_model, Char(46), :jl),
))

validate_free_tensors()
validate_simple_quadratic_tensor()

open(
    joinpath(@__DIR__, string(:success, Char(46), :sentinel));
    write=true,
    create=true,
    truncate=true,
) do io
    write(io, string(:ok))
end
