import Pkg

Pkg.activate(@__DIR__)
include(joinpath(
    @__DIR__,
    string(:examples),
    string(:Single_Flavor_Gross_Neveu_Wilson, Char(46), :jl),
))

validate_single_flavor_gnw_tensor()

open(
    joinpath(@__DIR__, string(:success, Char(46), :sentinel));
    write=true,
    create=true,
    truncate=true,
) do io
    write(io, string(:ok))
end
