##############################################################
##  GrassmannGenerator — a single anticommuting generator   ##
##############################################################

"""
    GrassmannGenerator

Represents a single Grassmann generator (an anticommuting variable).

Fields
- `name::Symbol` , e.g. `:η`, `:θ`, `:ψ`
- `index::Union{Nothing, Int}`, `nothing` means no index,
   positive integers give η₁, η₂, ...
- `bar::Bool`, `false` → η;  `true` → η̄  (conjugate / dual)

Two generators are considered equal iff all three fields match.
Any two distinct generators strictly anticommute.
"""

struct GrassmannGenerator
    name::Symbol
    index::Union{Nothing,Int}
    bar::Bool
end

# ##################### Constructors ############################

"""
    grassmann(name::Symbol; index=nothing, bar=false) -> GrassmannGenerator

Create a single Grassmann generator.

Example:
    η  = grassmann(:η)
    η̄  = grassmann(:η; bar=true)
    η1 = grassmann(:η; index=1)
"""

grassmann(name::Symbol; index=nothing, bar=false) = GrassmannGenerator(name, index, bar)

"""
    grassmann(name::Symbol, indices) -> Vector{GrassmannGenerator}

Create multiple indexed generators η₁, η₂, … from a range or collection.

Example:
    grassmann(:η, 1:3)  # returns [η₁, η₂, η₃]
"""

function grassmann(name::Symbol, inds_range)
    [GrassmannGenerator(name, i, false) for i in inds_range]
end

"""
    dual(g::GrassmannGenerator) -> GrassmannGenerator

Return the dual (bar-conjugate) of `g`. 

Example:
    dual(η) = η̄, `dual(η̄) = η`.
"""

dual(g::GrassmannGenerator) = GrassmannGenerator(g.name, g.index, !g.bar)

# ##################### Canonical ordering #####################

# We need a total order on generators so that every product can be written
# in a unique canonical (sorted) form.  The ordering groups bar/non-bar pairs
# together within the same index level:
#
#   1. name  (lexicographic)
#   2. index (nothing = unindexed sorts first, then by value)
#   3. bar   (false < true, so η before η̄ within the same index group)
#
# Example:  η < η̄ < η₁ < η₁̄ < η₂ < η₂̄ < θ < θ̄ < θ₁ < ψ  ...

function Base.isless(a::GrassmannGenerator, b::GrassmannGenerator)
    # 1. name
    sa, sb = String(a.name), String(b.name)
    sa != sb && return sa < sb
    # 2. index: treat nothing as -1 (sorts before any non-negative index)
    ai = a.index === nothing ? -1 : a.index
    bi = b.index === nothing ? -1 : b.index
    ai != bi && return ai < bi
    # 3. bar: false < true  (η before η̄ within the same (name, index) group)
    return !a.bar & b.bar
end

Base.isequal(a::GrassmannGenerator, b::GrassmannGenerator) = a.name == b.name && a.index == b.index && a.bar == b.bar

Base.:(==)(a::GrassmannGenerator, b::GrassmannGenerator) = isequal(a, b)

function Base.hash(g::GrassmannGenerator, h::UInt)
    h = hash(g.name, h)
    h = hash(g.index, h)
    h = hash(g.bar, h)
    return h
end

# ##################### @grassmann macro #####################

"""
    @grassmann η θ ψ ...

Convenience macro that creates one `GrassmannGenerator` for each symbol and
binds it to a local variable with the same name.

Example:
    @grassmann η θ ψ
    η * θ   
"""

macro grassmann(syms...)
    exprs = [:($(esc(s)) = GrassmannGenerator($(QuoteNode(s)), nothing, false))
             for s in syms]
    return Expr(:block, exprs...)
end
