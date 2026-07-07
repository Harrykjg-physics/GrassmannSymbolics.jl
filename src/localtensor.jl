##############################################################
##  Local Grassmann tensor construction                      ##
##############################################################

_as_grassmann_expr(e::GrassmannExpr) = e
_as_grassmann_expr(g::GrassmannGenerator) = GrassmannExpr(g)
_as_grassmann_expr(scalar) = GrassmannExpr(scalar)

"""
    integrate(integrand, measure_order) -> GrassmannExpr

Perform a Berezin integral with the measure written in `measure_order`.
For example, `integrate(f, [ψ₁, ψ̄₁, ψ₂, ψ̄₂])` represents
`∫dψ₁ dψ̄₁ dψ₂ dψ̄₂ f`; as usual, the rightmost differential acts first.

Keeping the measure order explicit is essential when comparing a generated
tensor with a convention fixed in the literature.
"""
function integrate(
    integrand,
    measure_order::AbstractVector{GrassmannGenerator},
)
    measure = BerezinMeasure(collect(measure_order))
    return measure * _as_grassmann_expr(integrand)
end

# Contract oriented pairs with d(barred)d(unbarred) exp(-barred*unbarred).
function contract_grassmann(tensors, pairs; simplify::Bool=true)
    integrand = GrassmannExpr(1)
    for tensor in tensors
        integrand = integrand * _as_grassmann_expr(tensor)
    end

    measure_order = GrassmannGenerator[]
    for (barred, unbarred) in pairs
        integrand = integrand * exp(-barred * unbarred)
        append!(measure_order, (barred, unbarred))
    end

    contracted = integrate(integrand, measure_order)
    return simplify ? simplify_coefficients(contracted) : contracted
end

"""
    local_grassmann_tensor(onsite_action, bond_exponents, measure_order;
                           simplify=true) -> GrassmannExpr

Construct the fundamental tensor at one lattice site,

    ∫ Dψ exp(-S_onsite) ∏ₖ exp(Bₖ),

where `bond_exponents` are the local factors produced by decomposing hopping
terms with auxiliary Grassmann fields. The original site fields are integrated
in `measure_order`; auxiliary fields remain as the tensor legs.

This is the common symbolic kernel for nearest-neighbor fermion models. A
model-specific compiler only has to turn its hopping matrices or interaction
terms into the local exponents `Bₖ`.
"""
function local_grassmann_tensor(
    onsite_action,
    bond_exponents,
    measure_order::AbstractVector{GrassmannGenerator};
    simplify::Bool=true,
)
    integrand = exp(-_as_grassmann_expr(onsite_action))
    for bond_exponent in bond_exponents
        integrand = integrand * exp(_as_grassmann_expr(bond_exponent))
    end
    tensor = integrate(integrand, measure_order)
    return simplify ? simplify_coefficients(tensor) : tensor
end
