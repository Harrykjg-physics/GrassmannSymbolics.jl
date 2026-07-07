##############################################################
##  GrassmannExpr — a general element of the exterior algebra ##
##############################################################

# A Monomial is a sorted (canonical) list of distinct generators.
# The empty vector represents the scalar "1".
const Monomial = Vector{GrassmannGenerator}

"""
    GrassmannExpr

An element of the exterior (Grassmann) algebra over an arbitrary coefficient
ring. Internally stored as a sparse dictionary

    terms :: Dict{Monomial, Any}

where each key is a *sorted*, duplicate-free `Vector{GrassmannGenerator}` (a
monomial) and the value is its coefficient (any Julia number or symbolic
expression from Symbolics.jl).

The scalar "1" corresponds to the empty monomial `GrassmannGenerator[]`.
"""
struct GrassmannExpr
    terms::Dict{Monomial,Any}
end

##############################################################
##  Internal helpers                                         ##
##############################################################

"""
    _iszero_coeff(c) -> Bool

Return `true` if `c` is numerically or symbolically zero.
"""
function _iszero_coeff(c)::Bool
    try
        return iszero(c)
    catch
        return false
    end
end

"""
    _sort_and_sign(gens::Vector{GrassmannGenerator}) -> (Monomial, Int)

Sort `gens` using insertion sort, counting the number of transpositions
needed.  Returns the sorted monomial and the sign `±1`.

If any generator appears more than once the result is `(Monomial[], 0)`,
signalling a zero contribution (nilpotency: g*g = 0).
"""
function _sort_and_sign(gens::Vector{GrassmannGenerator})
    n = length(gens)
    n == 0 && return (GrassmannGenerator[], 1)

    arr = copy(gens)
    sign = 1
    for i in 2:n
        j = i
        while j > 1
            if arr[j] < arr[j-1]
                arr[j], arr[j-1] = arr[j-1], arr[j]
                sign = -sign
                j -= 1
            elseif arr[j] == arr[j-1]
                # repeated generator → nilpotent, result is 0
                return (GrassmannGenerator[], 0)
            else
                break
            end
        end
    end
    return (arr, sign)
end

"""
    _cleanup!(expr::GrassmannExpr) -> GrassmannExpr

Remove all terms whose coefficient is zero (in-place).  Returns `expr`.
"""
function _cleanup!(expr::GrassmannExpr)
    for (k, v) in expr.terms
        if _iszero_coeff(v)
            delete!(expr.terms, k)
        end
    end
    return expr
end

##############################################################
##  Constructors                                             ##
##############################################################

# Internal: wrap a Dict directly (no copy)
_gexpr(d::Dict{Monomial,Any}) = GrassmannExpr(d)

"""
    GrassmannExpr(c::Number) -> GrassmannExpr

Wrap a scalar (including Symbolics.Num) as a GrassmannExpr.
"""

function GrassmannExpr(c::Number)
    _iszero_coeff(c) && return GrassmannExpr(Dict{Monomial,Any}())
    GrassmannExpr(Dict{Monomial,Any}(GrassmannGenerator[] => c))
end

"""
    GrassmannExpr(g::GrassmannGenerator) -> GrassmannExpr

Lift a single generator into the algebra.
"""

function GrassmannExpr(g::GrassmannGenerator)
    GrassmannExpr(Dict{Monomial,Any}(GrassmannGenerator[g] => 1))
end

##############################################################
##  Predicates                                               ##
##############################################################

Base.iszero(e::GrassmannExpr) = isempty(e.terms)

"""Return the scalar (grade-0) part of `e`."""
function scalar_part(e::GrassmannExpr)
    get(e.terms, GrassmannGenerator[], 0)
end

"""Return the pure-Grassmann part (all terms with grade ≥ 1)."""
function grassmann_part(e::GrassmannExpr)
    d = Dict{Monomial,Any}()
    for (k, v) in e.terms
        isempty(k) || (d[k] = v)
    end
    GrassmannExpr(d)
end

"""
    is_even(e::GrassmannExpr) -> Bool

Return `true` if every non-zero term has even grade.
"""

function is_even(e::GrassmannExpr)
    for (k, v) in e.terms
        _iszero_coeff(v) && continue
        iseven(length(k)) || return false
    end
    return true
end

"""
    is_odd(e::GrassmannExpr) -> Bool

Return `true` if every non-zero term has odd grade.
"""

function is_odd(e::GrassmannExpr)
    for (k, v) in e.terms
        _iszero_coeff(v) && continue
        isodd(length(k)) || return false
    end
    return true
end

"""
    is_grassmann(e::GrassmannExpr) -> Bool

Return `true` if `e` contains at least one term with grade ≥ 1.
"""

function is_grassmann(e::GrassmannExpr)
    for (k, v) in e.terms
        _iszero_coeff(v) && continue
        isempty(k) || return true
    end
    return false
end

"""
    generators(e::GrassmannExpr) -> Vector{GrassmannGenerator}

Return the sorted list of all generators appearing in `e`.
"""

# A single generator is a homogeneous odd Grassmann expression.
is_even(::GrassmannGenerator) = false
is_odd(::GrassmannGenerator) = true
is_grassmann(::GrassmannGenerator) = true

function generators(e::GrassmannExpr)
    unique_generators = Set{GrassmannGenerator}()
    for monomial in keys(e.terms)
        union!(unique_generators, monomial)
    end
    return sort!(collect(unique_generators))
end

##############################################################
##  Addition / Subtraction                                   ##
##############################################################

function Base.:+(a::GrassmannExpr, b::GrassmannExpr)
    d = copy(a.terms)
    for (k, v) in b.terms
        if haskey(d, k)
            d[k] = d[k] + v
        else
            d[k] = v
        end
    end
    _cleanup!(GrassmannExpr(d))
end

function Base.:+(a::GrassmannExpr, b)
    _iszero_coeff(b) && return a
    d = copy(a.terms)
    key = GrassmannGenerator[]
    d[key] = get(d, key, 0) + b
    _cleanup!(GrassmannExpr(d))
end

Base.:+(a, b::GrassmannExpr) = b + a
# Explicit overloads to resolve GrassmannGenerator vs GrassmannExpr ambiguity
Base.:+(g::GrassmannGenerator, b::GrassmannExpr) = GrassmannExpr(g) + b
Base.:+(a::GrassmannExpr, g::GrassmannGenerator) = a + GrassmannExpr(g)
Base.:+(g::GrassmannGenerator, b) = GrassmannExpr(g) + b
Base.:+(a, g::GrassmannGenerator) = a + GrassmannExpr(g)
Base.:+(g::GrassmannGenerator, h::GrassmannGenerator) =
    GrassmannExpr(g) + GrassmannExpr(h)

Base.:-(g::GrassmannGenerator) = -GrassmannExpr(g)

function Base.:-(e::GrassmannExpr)
    d = empty(e.terms)          # same Dict type, avoids spelling out Monomial
    for (k, v) in e.terms
        d[k] = -v
    end
    GrassmannExpr(d)
end
Base.:-(a::GrassmannExpr, b::GrassmannExpr) = a + (-b)
Base.:-(a::GrassmannExpr, b) = a + (-b)
Base.:-(a, b::GrassmannExpr) = a + (-b)
# Explicit overloads to resolve GrassmannGenerator vs GrassmannExpr ambiguity
Base.:-(g::GrassmannGenerator, b::GrassmannExpr) = GrassmannExpr(g) - b
Base.:-(a::GrassmannExpr, g::GrassmannGenerator) = a - GrassmannExpr(g)
Base.:-(a::GrassmannGenerator, b) = GrassmannExpr(a) - b
Base.:-(a, b::GrassmannGenerator) = a - GrassmannExpr(b)
Base.:-(a::GrassmannGenerator, b::GrassmannGenerator) =
    GrassmannExpr(a) - GrassmannExpr(b)

##############################################################
##  Multiplication                                           ##
##############################################################

"""
    _monomial_product(m1::Monomial, m2::Monomial) -> (Monomial, Int)

Compute the exterior (wedge) product of two sorted monomials.
Returns the resulting sorted monomial and the sign.
A sign of 0 means the result is zero (repeated generator).
"""

function _monomial_product(m1::Monomial, m2::Monomial)
    isempty(m1) && return (copy(m2), 1)
    isempty(m2) && return (copy(m1), 1)

    n1, n2 = length(m1), length(m2)
    merged = GrassmannGenerator[]
    sizehint!(merged, n1 + n2)
    i = j = 1
    odd_permutation = false

    while i <= n1 && j <= n2
        if m1[i] < m2[j]
            push!(merged, m1[i])
            i += 1
        elseif m2[j] < m1[i]
            push!(merged, m2[j])
            # Moving m2[j] before every unconsumed generator of m1 costs
            # n1-i+1 transpositions. Only its parity is needed.
            isodd(n1 - i + 1) && (odd_permutation = !odd_permutation)
            j += 1
        else
            return (GrassmannGenerator[], 0)
        end
    end

    while i <= n1
        push!(merged, m1[i])
        i += 1
    end
    while j <= n2
        push!(merged, m2[j])
        j += 1
    end

    return (merged, odd_permutation ? -1 : 1)
end

function Base.:*(a::GrassmannExpr, b::GrassmannExpr)
    d = Dict{Monomial,Any}()
    for (k1, v1) in a.terms
        for (k2, v2) in b.terms
            (mono, sgn) = _monomial_product(k1, k2)
            sgn == 0 && continue
            coeff = v1 * v2 * sgn
            if haskey(d, mono)
                d[mono] = d[mono] + coeff
            else
                d[mono] = coeff
            end
        end
    end
    _cleanup!(GrassmannExpr(d))
end

# scalar * GrassmannExpr
function Base.:*(a, b::GrassmannExpr)
    _iszero_coeff(a) && return GrassmannExpr(Dict{Monomial,Any}())
    GrassmannExpr(Dict{Monomial,Any}(k => a * v for (k, v) in b.terms))
end

Base.:*(a::GrassmannExpr, b) = b * a

# GrassmannGenerator participates in the algebra as a degree-1 element
Base.:*(g::GrassmannGenerator, h::GrassmannGenerator) =
    GrassmannExpr(g) * GrassmannExpr(h)
Base.:*(g::GrassmannGenerator, b::GrassmannExpr) = GrassmannExpr(g) * b
Base.:*(a::GrassmannExpr, g::GrassmannGenerator) = a * GrassmannExpr(g)
Base.:*(a, g::GrassmannGenerator) = a * GrassmannExpr(g)
Base.:*(g::GrassmannGenerator, b) = GrassmannExpr(g) * b

##############################################################
##  Division by scalar                                       ##
##############################################################

Base.:/(e::GrassmannExpr, s) =
    GrassmannExpr(Dict{Monomial,Any}(k => v / s for (k, v) in e.terms))

Base.:/(g::GrassmannGenerator, s) = GrassmannExpr(g) / s

##############################################################
##  Power                                                    ##
##############################################################

"""
    ^(e::GrassmannExpr, n::Integer) -> GrassmannExpr

Compute the n-th exterior power of `e`.  For odd-grade elements `e^2 = 0`;
the Taylor expansion terminates automatically.
"""
function Base.:^(e::GrassmannExpr, n::Integer)
    n < 0  && error("Negative powers of Grassmann expressions are not defined.")
    n == 0 && return GrassmannExpr(1)
    n == 1 && return e
    result = GrassmannExpr(1)
    base   = e
    k = n
    while k > 0
        if isodd(k)
            result = result * base
            # Short-circuit: once result is zero it stays zero
            iszero(result) && return result
        end
        k >>= 1
        base = base * base
    end
    return result
end

Base.:^(g::GrassmannGenerator, n::Integer) = GrassmannExpr(g)^n

##############################################################
##  exp                                                      ##
##############################################################

"""
    exp(e::GrassmannExpr) -> GrassmannExpr

Compute `exp(e)` using the finite Taylor expansion.

Since the pure-Grassmann part `X = e - scalar_part(e)` is nilpotent
(all generators satisfy g² = 0), the series terminates:

    exp(c + X) = exp(c) · (1 + X + X²/2! + … + Xⁿ/n!)

where `n` is the number of distinct generators in `e`.
The scalar `exp(c)` is computed by Julia dispatch (handles both
`Number` and `Symbolics.Num` automatically).
"""

function Base.exp(e::GrassmannExpr)
    c = scalar_part(e)
    X = grassmann_part(e)

    ngens = length(generators(e))
    ngens == 0 && return GrassmannExpr(exp(c))  # purely scalar

    # Taylor series: sum_{p=0}^{ngens} X^p / p!
    result  = GrassmannExpr(1)   # p = 0 term
    Xp      = GrassmannExpr(1)   # X^0
    pfact   = 1                  # 0!
    for p in 1:ngens
        pfact *= p
        Xp = Xp * X
        iszero(Xp) && break
        result = result + Xp / pfact
    end

    # multiply by exp(c)
    ec = exp(c)
    return ec * result
end

Base.exp(g::GrassmannGenerator) = exp(GrassmannExpr(g))

##############################################################
##  get_coeff                                                ##
##############################################################

"""
    get_coeff(e::GrassmannExpr, gens::Vector{GrassmannGenerator}) -> Any

Extract the coefficient of the monomial formed by `gens` in `e`.
`gens` need not be sorted; the canonical sign is applied automatically.

Returns 0 if that monomial does not appear.

# Example
```julia
@grassmann η θ
e = 3η*θ + 2η
get_coeff(e, [η, θ])   # → 3
get_coeff(e, [θ, η])   # → -3  (anticommutation sign)
```
"""

function get_coeff(e::GrassmannExpr, gens::Vector{GrassmannGenerator})
    isempty(gens) && return scalar_part(e)
    (mono, sgn) = _sort_and_sign(copy(gens))
    sgn == 0 && return 0
    coeff = get(e.terms, mono, 0)
    return sgn * coeff
end

##############################################################
##  Coefficient transformations                              ##
##############################################################

"""
    map_coefficients(f, e::GrassmannExpr) -> GrassmannExpr

Apply `f` independently to every scalar coefficient of `e`. Grassmann
monomials and their canonical ordering are left unchanged.
"""
function map_coefficients(f, e::GrassmannExpr)
    result = Dict{Monomial,Any}()
    for (monomial, coefficient) in e.terms
        mapped = f(coefficient)
        _iszero_coeff(mapped) || (result[copy(monomial)] = mapped)
    end
    return GrassmannExpr(result)
end

"""
    simplify_coefficients(e::GrassmannExpr; expand=true) -> GrassmannExpr

Simplify every scalar coefficient with Symbolics.jl. This is useful after a
local Grassmann integration, where algebraically identical tensor entries can
have different unsimplified expression trees.
"""
function simplify_coefficients(e::GrassmannExpr; expand::Bool=true)
    simplify_one(coefficient) = try
        Symbolics.simplify(coefficient; expand=expand)
    catch
        try
            Symbolics.simplify(coefficient)
        catch
            coefficient
        end
    end
    return map_coefficients(simplify_one, e)
end

"""
    substitute_coefficients(e::GrassmannExpr, rules) -> GrassmannExpr

Apply `Symbolics.substitute` to every coefficient. Numeric coefficients and
other coefficient types unsupported by Symbolics are kept unchanged.
"""
function substitute_coefficients(e::GrassmannExpr, rules)
    substitute_one(coefficient) = try
        Symbolics.substitute(coefficient, rules)
    catch
        coefficient
    end
    return map_coefficients(substitute_one, e)
end

# Relabel generators and restore canonical order, including its sign.
function map_generators(f, e::GrassmannExpr)
    result = Dict{Monomial,Any}()
    for (monomial, coefficient) in e.terms
        mapped = GrassmannGenerator[f(generator) for generator in monomial]
        canonical, sign = _sort_and_sign(mapped)
        sign == 0 && continue
        result[canonical] = get(result, canonical, 0) + sign * coefficient
    end
    return _cleanup!(GrassmannExpr(result))
end

# Rename generators in rules while preserving every absent generator.
function substitute_generators(e::GrassmannExpr, rules)
    return map_generators(generator -> get(rules, generator, generator), e)
end
"""
    symbolically_equal(a::GrassmannExpr, b::GrassmannExpr) -> Bool

Return `true` when all coefficients of `a-b` simplify to zero.
"""
function symbolically_equal(a::GrassmannExpr, b::GrassmannExpr)
    difference = simplify_coefficients(a - b)
    return all(_iszero_coeff, values(difference.terms))
end
##############################################################
##  Equality                                                 ##
##############################################################

function Base.:(==)(a::GrassmannExpr, b::GrassmannExpr)
    d = a - b
    return iszero(d)
end
