using GrassmannSymbolics
using Symbolics
include(joinpath(@__DIR__, "..", "examples", "Free_Wilson_and_Staggered.jl"))

function rephase(expression, phases)
    result = GrassmannExpr(0)
    for (monomial, coefficient) in expression.terms
        factor = prod((get(phases, generator.name, 1) for generator in monomial); init=1)
        term = GrassmannExpr(factor * coefficient)
        for generator in monomial
            term = term * generator
        end
        result = result + term
    end
    return result
end

function substitute_coefficients(expression, rules)
    return map_coefficients(expression) do coefficient
        try
            Symbolics.substitute(coefficient, rules)
        catch
            coefficient
        end
    end
end

function numeric_value(coefficient)
    value = try
        Symbolics.value(coefficient)
    catch
        coefficient
    end
    return Float64(value)
end

function max_coefficient_error(left, right, rules)
    difference = simplify_coefficients(
        substitute_coefficients(left - right, rules),
    )
    isempty(difference.terms) && return 0.0
    return maximum(abs(numeric_value(c)) for c in values(difference.terms))
end

@variables D A11 A12 A21 A22 B11 B22 M p
q = inv(sqrt(big"2"))
wilson_rules = Dict(
    D => big"1.375",
    A11 => -q,
    A12 => -q,
    A21 => -q,
    A22 => q,
    B11 => big"-1",
    B22 => big"-1",
)

println("Wilson phase/channel search")
for channels in ((1, 2), (2, 1))
    result = derive_wilson_tensor(
        x_channel=channels[1],
        y_channel=channels[2],
    )
    best = (Inf, nothing)
    for signs in Iterators.product((-1, 1), (-1, 1), (-1, 1), (-1, 1))
        phases = Dict(
            :ηx => signs[1],
            :ζx => signs[2],
            :ηt => signs[3],
            :ζt => signs[4],
        )
        error = max_coefficient_error(
            rephase(result.tensor, phases),
            result.published_tensor,
            wilson_rules,
        )
        error < best[1] && (best = (error, copy(phases)))
    end
    println("  channels=", channels, " best_error=", best[1], " phases=", best[2])
end

println("Staggered phase search")
staggered = derive_staggered_tensor()
for parity in (-1, 1)
    rules = Dict(M => big"0.625", p => parity)
    best = (Inf, nothing)
    for signs in Iterators.product((-1, 1), (-1, 1), (-1, 1), (-1, 1))
        phases = Dict(
            :ηx => signs[1],
            :ζx => signs[2],
            :ηt => signs[3],
            :ζt => signs[4],
        )
        error = max_coefficient_error(
            rephase(staggered.tensor, phases),
            staggered.published_tensor,
            rules,
        )
        error < best[1] && (best = (error, copy(phases)))
    end
    println("  p=", parity, " best_error=", best[1], " phases=", best[2])
end

println("Staggered one-site closure at p=1")
derived_p1 = specialize_staggered_parity(staggered.tensor, p, 1)
published_p1 = specialize_staggered_parity(staggered.published_tensor, p, 1)
println("  action-derived closure = ", contract_periodic_one_site(derived_p1, staggered.legs))
println("  published closure      = ", contract_periodic_one_site(published_p1, staggered.legs))
println("  direct onsite result   = -M (hopping cancels on a self-link)")
# From (B.6), the minor with rows (ηx, ζbar_t) and columns
# (ηbar_x, ηbar_t) equals 1/2, so its bilinear coefficient matrix
# cannot be the rank-one outer product implied by the one-component
# local integral in equation (3.10).
minor_B6 = (-1 // 2) * (-1 // 2) - (-1 // 2) * (1 // 2)
println("B.6 diagnostic 2x2 minor = ", minor_B6)