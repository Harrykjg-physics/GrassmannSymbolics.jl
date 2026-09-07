##############################################################
##  Local Grassmann tensor construction                      ##
##############################################################

_as_grassmann_expr(e::GrassmannExpr) = e
_as_grassmann_expr(g::GrassmannGenerator) = GrassmannExpr(g)
_as_grassmann_expr(scalar) = GrassmannExpr(scalar)

struct LocalChannel
    exponent::Any
    legs::Vector{GrassmannGenerator}
    tag::Union{Nothing,Symbol}
end

function LocalChannel(
    exponent;
    legs::AbstractVector{GrassmannGenerator}=GrassmannGenerator[],
    tag::Union{Nothing,Symbol}=nothing,
)
    return LocalChannel(exponent, collect(legs), tag)
end

function LocalChannel(
    exponent,
    leg::GrassmannGenerator;
    tag::Union{Nothing,Symbol}=nothing,
)
    return LocalChannel(exponent, [leg], tag)
end

struct NearestNeighborTensorSpec
    onsite_action::Any
    measure_order::Vector{GrassmannGenerator}
    channels::Vector{LocalChannel}
    leg_order::Vector{GrassmannGenerator}
    metadata::Dict{Symbol,Any}
end

function _infer_leg_order(channels)
    legs = GrassmannGenerator[]
    for channel in channels
        for leg in channel.legs
            leg in legs || push!(legs, leg)
        end
    end
    return legs
end

function NearestNeighborTensorSpec(
    onsite_action,
    measure_order::AbstractVector{GrassmannGenerator},
    channels::AbstractVector{LocalChannel};
    leg_order::Union{Nothing,AbstractVector{GrassmannGenerator}}=nothing,
    metadata=Dict{Symbol,Any}(),
)
    ordered_legs = leg_order === nothing ? _infer_leg_order(channels) : collect(leg_order)
    return NearestNeighborTensorSpec(
        onsite_action,
        collect(measure_order),
        collect(channels),
        ordered_legs,
        Dict{Symbol,Any}(metadata),
    )
end

function split_hopping_channel(
    left_factor,
    auxiliary::GrassmannGenerator,
    right_factor;
    left_scale=1,
    right_scale=1,
    tag::Union{Nothing,Symbol}=nothing,
)
    auxiliary_bar = dual(auxiliary)
    return LocalChannel[
        LocalChannel(left_scale * left_factor * auxiliary, auxiliary; tag=tag),
        LocalChannel(right_scale * auxiliary_bar * right_factor, auxiliary_bar; tag=tag),
    ]
end


abstract type LocalActionTerm end

struct OnsiteTerm <: LocalActionTerm
    action::Any
    tag::Union{Nothing,Symbol}
end

OnsiteTerm(action; tag::Union{Nothing,Symbol}=nothing) = OnsiteTerm(action, tag)

struct FactorizedHoppingTerm <: LocalActionTerm
    left_factor::Any
    auxiliary::GrassmannGenerator
    right_factor::Any
    left_scale::Any
    right_scale::Any
    tag::Union{Nothing,Symbol}
end

function FactorizedHoppingTerm(
    left_factor,
    auxiliary::GrassmannGenerator,
    right_factor;
    left_scale=1,
    right_scale=1,
    tag::Union{Nothing,Symbol}=nothing,
)
    return FactorizedHoppingTerm(
        left_factor,
        auxiliary,
        right_factor,
        left_scale,
        right_scale,
        tag,
    )
end

struct NearestNeighborAction
    measure_order::Vector{GrassmannGenerator}
    terms::Vector{LocalActionTerm}
    channels::Vector{LocalChannel}
    leg_order::Union{Nothing,Vector{GrassmannGenerator}}
    metadata::Dict{Symbol,Any}
end

function NearestNeighborAction(
    measure_order::AbstractVector{GrassmannGenerator};
    terms::AbstractVector{<:LocalActionTerm}=LocalActionTerm[],
    onsite_terms::AbstractVector=Any[],
    hopping_terms::AbstractVector{<:FactorizedHoppingTerm}=FactorizedHoppingTerm[],
    channels::AbstractVector{LocalChannel}=LocalChannel[],
    leg_order::Union{Nothing,AbstractVector{GrassmannGenerator}}=nothing,
    metadata=Dict{Symbol,Any}(),
)
    all_terms = LocalActionTerm[]
    append!(all_terms, terms)
    append!(all_terms, OnsiteTerm.(onsite_terms))
    append!(all_terms, hopping_terms)
    ordered_legs = leg_order === nothing ? nothing : collect(leg_order)
    return NearestNeighborAction(
        collect(measure_order),
        all_terms,
        collect(channels),
        ordered_legs,
        Dict{Symbol,Any}(metadata),
    )
end

function nearest_neighbor_tensor_spec(action::NearestNeighborAction)
    onsite_action = GrassmannExpr(0)
    channels = copy(action.channels)
    term_tags = Symbol[]

    for term in action.terms
        if term isa OnsiteTerm
            onsite_action = onsite_action + _as_grassmann_expr(term.action)
            term.tag === nothing || push!(term_tags, term.tag)
        elseif term isa FactorizedHoppingTerm
            append!(
                channels,
                split_hopping_channel(
                    term.left_factor,
                    term.auxiliary,
                    term.right_factor;
                    left_scale=term.left_scale,
                    right_scale=term.right_scale,
                    tag=term.tag,
                ),
            )
            term.tag === nothing || push!(term_tags, term.tag)
        else
            throw(ArgumentError("unsupported local action term"))
        end
    end

    metadata = copy(action.metadata)
    isempty(term_tags) || (metadata[:term_tags] = term_tags)
    return NearestNeighborTensorSpec(
        onsite_action,
        action.measure_order,
        channels;
        leg_order=action.leg_order,
        metadata=metadata,
    )
end

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

function compile_nearest_neighbor_tensor(
    spec::NearestNeighborTensorSpec;
    simplify::Bool=true,
)
    tensor = local_grassmann_tensor(
        spec.onsite_action,
        [channel.exponent for channel in spec.channels],
        spec.measure_order;
        simplify=simplify,
    )
    return (;
        tensor,
        legs=copy(spec.leg_order),
        channels=copy(spec.channels),
        measure_order=copy(spec.measure_order),
        metadata=copy(spec.metadata),
    )
end

function compile_nearest_neighbor_tensor(
    action::NearestNeighborAction;
    simplify::Bool=true,
)
    spec = nearest_neighbor_tensor_spec(action)
    return compile_nearest_neighbor_tensor(spec; simplify=simplify)
end

function compile_nearest_neighbor_tensor(
    onsite_action,
    measure_order::AbstractVector{GrassmannGenerator},
    channels::AbstractVector{LocalChannel};
    leg_order=nothing,
    metadata=Dict{Symbol,Any}(),
    simplify::Bool=true,
)
    spec = NearestNeighborTensorSpec(
        onsite_action,
        measure_order,
        channels;
        leg_order=leg_order,
        metadata=metadata,
    )
    return compile_nearest_neighbor_tensor(spec; simplify=simplify)
end
