##############################################################
##  Unicode display for GrassmannGenerator, GrassmannExpr   ##
##  and BerezinMeasure                                       ##
##############################################################

# Map digit characters to Unicode subscript equivalents
const _SUBSCRIPT_DIGITS = Dict(
    '0' => '₀', '1' => '₁', '2' => '₂', '3' => '₃',
    '4' => '₄', '5' => '₅', '6' => '₆', '7' => '₇',
    '8' => '₈', '9' => '₉',
)

"""Turn an integer into a Unicode subscript string, e.g. 12 → "₁₂"."""
function _subscript(n::Int)
    s = string(n)
    join(get(_SUBSCRIPT_DIGITS, c, c) for c in s)
end

##############################################################
##  GrassmannGenerator                                       ##
##############################################################

function Base.show(io::IO, g::GrassmannGenerator)
    print(io, _generator_label(g))
end

"""Return the display string for a single generator."""
function _generator_label(g::GrassmannGenerator)
    base = String(g.name)
    idx  = g.index === nothing ? "" : _subscript(g.index)
    bar  = g.bar ? "̄" : ""   # combining overbar U+0304
    return base * idx * bar
end

##############################################################
##  GrassmannExpr                                            ##
##############################################################

"""
Return the display string for a sorted monomial (product of generators).
Empty monomial → "1".
"""
function _monomial_label(mono::Monomial)
    isempty(mono) && return "1"
    join(_generator_label(g) for g in mono)
end

function Base.show(io::IO, e::GrassmannExpr)
    if isempty(e.terms)
        print(io, "0")
        return
    end

    # Sort terms for deterministic output:
    # first by grade, then lexicographically by monomial label
    sorted_pairs = sort(
        collect(e.terms),
        by = p -> (length(p.first), _monomial_label(p.first)),
    )

    parts = String[]
    for (mono, coeff) in sorted_pairs
        _iszero_coeff(coeff) && continue
        mlabel = _monomial_label(mono)
        clabel = _coeff_label(coeff, mlabel)
        push!(parts, clabel)
    end

    isempty(parts) && (print(io, "0"); return)

    # Join: first term as-is, remaining prefixed with " + " / " - "
    buf = IOBuffer()
    print(buf, parts[1])
    for p in parts[2:end]
        if startswith(p, "-")
            print(buf, " - ", p[nextind(p, 1):end])
        else
            print(buf, " + ", p)
        end
    end
    print(io, String(take!(buf)))
end

"""Format a coefficient together with the monomial label it multiplies."""
function _coeff_label(coeff, mlabel::String)
    if mlabel == "1"
        return string(coeff)
    end
    # Try to detect simple ±1 coefficients numerically
    try
        if coeff == 1
            return mlabel
        elseif coeff == -1
            return "-" * mlabel
        end
    catch
    end
    return "(" * string(coeff) * ")" * mlabel
end

function Base.show(io::IO, ::MIME"text/plain", e::GrassmannExpr)
    show(io, e)
end

##############################################################
##  BerezinMeasure                                           ##
##############################################################

function Base.show(io::IO, m::BerezinMeasure)
    if isempty(m.generators)
        print(io, "∫(trivial)")
        return
    end
    for g in m.generators
        print(io, "∫d", _generator_label(g), " ")
    end
end

function Base.show(io::IO, ::MIME"text/plain", m::BerezinMeasure)
    show(io, m)
end
