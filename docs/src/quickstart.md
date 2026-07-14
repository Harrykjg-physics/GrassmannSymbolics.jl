# Quick start

Load the package and create generators:

```julia
julia> using GrassmannSymbolics

julia> @grassmann η ξ

julia> η * ξ
ηξ

julia> ξ * η
-ηξ

julia> η * η
0
```

Create indexed generators and their duals:

```julia
ψ = grassmann(:ψ, 1:2)
ψbar = dual.(ψ)
```

Finite Grassmann exponentials terminate automatically:

```julia
julia> exp(η * ξ)
1 + ηξ
```

Berezin integration uses explicit measures:

```julia
julia> integrate(η * ξ, [η, ξ])
-1

julia> integrate(η * ξ, [ξ, η])
1
```

The two examples differ only by the written measure order.  This is deliberate:
published Grassmann tensor formulas often differ by precisely these sign
conventions.

## Symbolic coefficients

`GrassmannSymbolics.jl` works with `Symbolics.jl` coefficients:

```julia
using Symbolics
@variables m g
@grassmann ψ
ψbar = dual(ψ)

expr = exp(m * ψbar * ψ)
# 1 + m * ψbar * ψ, up to the package's canonical generator order
```

Coefficient transformations are coefficient-wise:

```julia
simplify_coefficients(expr)
substitute_coefficients(expr, Dict(m => 2.0))
```

## A one-site tensor

```julia
@grassmann χ η
χbar = dual(χ)
ηbar = dual(η)

T = local_grassmann_tensor(
    3χbar * χ,
    [-χbar * η + (1 // 2) * ηbar * χ],
    [χ, χbar],
)
```

The physical fields `χ, χbar` are integrated out.  The auxiliary fields
`η, ηbar` remain as tensor legs.
