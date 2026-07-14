# Grassmann algebra

## Generators

The basic object is `GrassmannGenerator`, created with:

```julia
η = grassmann(:η)
ψ = grassmann(:ψ, 1:4)
ηbar = dual(η)
```

Each generator has three pieces of identity:

- `name::Symbol`
- `index::Union{Nothing,Int}`
- `bar::Bool`

Two distinct generators anticommute.  Repeating a generator in a product gives
zero.

## Expressions

`GrassmannExpr` stores a sparse linear combination of canonical monomials.  The
empty monomial represents the scalar part.  Coefficients may be ordinary Julia
numbers or symbolic objects such as `Symbolics.Num`.

Common operations:

```julia
scalar_part(expr)
grassmann_part(expr)
generators(expr)
is_even(expr)
is_odd(expr)
is_grassmann(expr)
```

## Coefficient extraction

Use `get_coeff(expr, gens)` to read the coefficient of a monomial.  The order in
`gens` matters: if it differs from canonical order, the anticommutation sign is
included.

```julia
@grassmann η ξ
expr = 3η * ξ

get_coeff(expr, [η, ξ])  #  3
get_coeff(expr, [ξ, η])  # -3
```

## Berezin integration

`integrate(expr, measure_order)` applies the Berezin measure in the written
order.  The rightmost differential acts first, so changing `measure_order` may
change a sign.

```julia
integrate(η * ξ, [ξ, η])  # 1
integrate(η * ξ, [η, ξ])  # -1
```

Keeping this convention explicit is essential for comparing against tensor
network references.
