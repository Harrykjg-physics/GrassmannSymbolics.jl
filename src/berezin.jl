##############################################################
##  BerezinMeasure — Berezin integral operator               ##
##############################################################

"""
    BerezinMeasure

Represents an ordered Berezin integration measure

    ∫ dg₁ dg₂ … dgₖ

The generators are stored in the order they appear *from left to right*
in the integral sign, which is the order in which integration is applied
(innermost first, i.e. `generators[end]` is applied first).

# Construction
Use [`d`](@ref) to build measures:

```julia
@grassmann η θ
d(η)          # ∫dη
d(η) * d(θ)   # ∫dη dθ  (apply dθ first, then dη)
```
"""
struct BerezinMeasure
    generators::Vector{GrassmannGenerator}
end

##############################################################
##  Constructor helper                                       ##
##############################################################

"""
    d(g::GrassmannGenerator) -> BerezinMeasure

Create the elementary Berezin measure `∫dg`.

Combine measures with `*`:

```julia
d(η) * d(θ)    # ∫dη dθ
```
"""
d(g::GrassmannGenerator) = BerezinMeasure([g])

##############################################################
##  Combining measures                                       ##
##############################################################

Base.:*(a::BerezinMeasure, b::BerezinMeasure) =
    BerezinMeasure(vcat(a.generators, b.generators))

##############################################################
##  Applying the measure to a GrassmannExpr                  ##
##############################################################

"""
    *(measure::BerezinMeasure, integrand) -> GrassmannExpr or scalar

Apply the Berezin integral to `integrand`.

### Rules

- Generators are applied **right to left** (innermost first):
  `d(η) * d(θ) * f` applies `∫dθ` first, then `∫dη`.
- The elementary rule is:  `∫dg (g) = 1`,  `∫dg (1) = 0`.
- If a generator appears more than once in the measure, the result is 0.

### Examples

```julia
@grassmann η θ ψ

d(η) * (η)            # → 1
d(η) * (θ)            # → 0
d(η) * (η * θ)        # → θ
d(θ) * d(η) * (η * θ) # → 1    (dη applied first, then dθ)
d(η) * d(θ) * (η * θ) # → −1   (dθ applied first, sign from ordering)
```
"""
function Base.:*(measure::BerezinMeasure, integrand::GrassmannExpr)
    # Check for repeated generators in the measure itself → result is 0
    mgens = measure.generators
    if length(unique(mgens)) < length(mgens)
        return GrassmannExpr(Dict{Monomial,Any}())
    end

    # Apply integrations right to left
    result = integrand
    for g in reverse(mgens)
        result = _apply_single_integral(g, result)
        iszero(result) && return result
    end
    return result
end

Base.:*(measure::BerezinMeasure, integrand::GrassmannGenerator) =
    measure * GrassmannExpr(integrand)

Base.:*(measure::BerezinMeasure, s) =
    _iszero_coeff(s) ? GrassmannExpr(Dict{Monomial,Any}()) :
    error("BerezinMeasure can only be applied to a GrassmannExpr or scalar 0.")

##############################################################
##  Single-generator Berezin integral (internal)             ##
##############################################################

"""
    _apply_single_integral(g::GrassmannGenerator, e::GrassmannExpr)
        -> GrassmannExpr

Apply the Berezin integral ∫dg to `e`.

For each term `c · g₁g₂…gₖ` in `e`:
- If `g` is not among the generators, the term vanishes.
- If `g` appears at position `i` (1-indexed in the sorted monomial),
  the sign is `(-1)^(i-1)` (from anti-commuting `g` to the leftmost position
  so that it can be cancelled by the measure), and the remaining monomial
  is `g₁…gᵢ₋₁ gᵢ₊₁…gₖ`.
"""
function _apply_single_integral(g::GrassmannGenerator, e::GrassmannExpr)
    d = Dict{Monomial,Any}()
    for (mono, coeff) in e.terms
        _iszero_coeff(coeff) && continue
        idx = findfirst(==(g), mono)
        idx === nothing && continue   # g not in this term → drops out

        # sign: moving g from position idx to the front costs (idx-1) swaps
        sgn = iseven(idx - 1) ? 1 : -1
        new_mono = deleteat!(copy(mono), idx)
        new_coeff = sgn * coeff
        if haskey(d, new_mono)
            d[new_mono] = d[new_mono] + new_coeff
        else
            d[new_mono] = new_coeff
        end
    end
    _cleanup!(GrassmannExpr(d))
end
