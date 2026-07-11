# Coefficient-array extraction utilities for Grassmann tensors.
#
# A Grassmann tensor is stored as a sparse GrassmannExpr.  Given an ordered
# grouping of external Grassmann legs, these routines materialize the scalar
# coefficient tensor as an ordinary Julia Array.  Each group becomes one Array
# dimension of size 2^length(group).  By default, grouped bits use the
# parity-preserving order common in Grassmann tensor-network codes: all even
# parity states first, then all odd parity states.

function _as_generator_vector(group)
    return GrassmannGenerator[generator for generator in group]
end

function _normalize_leg_groups(leg_groups)
    groups = collect(leg_groups)
    if isempty(groups)
        return Vector{GrassmannGenerator}[]
    end
    if all(group -> group isa GrassmannGenerator, groups)
        return [GrassmannGenerator[group] for group in groups]
    end
    return [_as_generator_vector(group) for group in groups]
end

function _validate_leg_groups(groups)
    flat = GrassmannGenerator[]
    for group in groups
        append!(flat, group)
    end
    length(unique(flat)) == length(flat) ||
        throw(ArgumentError("coefficient_array leg groups contain duplicate Grassmann generators"))
    return flat
end

function _leg_group_widths(leg_groups)
    groups = collect(leg_groups)
    if all(group -> group isa Integer, groups)
        return Int[group for group in groups]
    end
    return Int[length(group) for group in _normalize_leg_groups(leg_groups)]
end

coefficient_array_shape(leg_groups) = Tuple(2^width for width in _leg_group_widths(leg_groups))

function _lexicographic_occupation_bits(mask::Integer, width::Integer)
    return ntuple(k -> (mask >> (width - k)) & 1, width)
end

function occupation_bit_tuples(width::Integer; index_order::Symbol=:parity)
    width >= 0 || throw(ArgumentError("occupation bit width must be non-negative"))
    tuples = [_lexicographic_occupation_bits(mask, width) for mask in 0:(2^width - 1)]
    index_order === :parity && return vcat(
        [bits for bits in tuples if iseven(sum(bits))],
        [bits for bits in tuples if isodd(sum(bits))],
    )
    index_order === :lexicographic && return tuples
    index_order === :little_endian && return [
        ntuple(k -> (mask >> (k - 1)) & 1, width) for mask in 0:(2^width - 1)
    ]
    throw(ArgumentError("unknown coefficient_array index_order: $(index_order)"))
end

function occupation_bits_from_index(
    index::Integer,
    width::Integer;
    index_order::Symbol=:parity,
)
    index >= 1 || throw(ArgumentError("coefficient_array indices are 1-based"))
    ordered_bits = occupation_bit_tuples(width; index_order)
    index <= length(ordered_bits) ||
        throw(ArgumentError("index $index is too large for a $width-bit leg group"))
    return ordered_bits[index]
end

function occupation_index_from_bits(bits; index_order::Symbol=:parity)
    ordered_bits = occupation_bit_tuples(length(bits); index_order)
    index = findfirst(==(Tuple(bits)), ordered_bits)
    index === nothing && throw(ArgumentError("invalid occupation bits: $(bits)"))
    return index
end

function occupation_bits_from_indices(indices, leg_groups; index_order::Symbol=:parity)
    widths = _leg_group_widths(leg_groups)
    length(indices) == length(widths) ||
        throw(ArgumentError("number of indices must match number of leg groups"))
    chunks = [
        occupation_bits_from_index(index, width; index_order)
        for (index, width) in zip(indices, widths)
    ]
    isempty(chunks) && return ()
    return Tuple(vcat((collect(chunk) for chunk in chunks)...))
end

function occupied_legs_from_bits(leg_groups, bits)
    groups = _normalize_leg_groups(leg_groups)
    flat = _validate_leg_groups(groups)
    length(bits) == length(flat) ||
        throw(ArgumentError("number of occupation bits must match total number of Grassmann legs"))
    return GrassmannGenerator[leg for (leg, bit) in zip(flat, bits) if bit == 1]
end

function _coefficient_substitute(coefficient, substitutions)
    substitutions === nothing && return coefficient
    return try
        Symbolics.substitute(coefficient, substitutions)
    catch
        coefficient
    end
end

function _coefficient_simplify(coefficient, simplify::Bool)
    simplify || return coefficient
    return try
        Symbolics.simplify(coefficient; expand=true)
    catch
        try
            Symbolics.simplify(coefficient)
        catch
            coefficient
        end
    end
end

function _scalar_coefficient(value; substitutions=nothing, simplify::Bool=false)
    scalar = value isa GrassmannExpr ? scalar_part(value) : value
    scalar = _coefficient_substitute(scalar, substitutions)
    return _coefficient_simplify(scalar, simplify)
end

function _prepare_tensor_expression(tensor::GrassmannExpr, substitutions)
    substitutions === nothing && return tensor
    return substitute_coefficients(tensor, substitutions)
end

function _fill_coefficient_array!(
    array,
    coefficient_at_bits,
    leg_groups;
    substitutions=nothing,
    simplify::Bool=false,
    index_order::Symbol=:parity,
)
    for cartesian_index in CartesianIndices(array)
        bits = occupation_bits_from_indices(Tuple(cartesian_index), leg_groups; index_order)
        array[cartesian_index] = _scalar_coefficient(
            coefficient_at_bits(bits);
            substitutions=substitutions,
            simplify=simplify,
        )
    end
    return array
end

"""
    coefficient_array(tensor::GrassmannExpr, leg_groups; substitutions=nothing,
                      simplify=false, element_type=Any, index_order=:parity)

Materialize the scalar coefficient tensor of a Grassmann expression as a Julia
`Array`.

`leg_groups` specifies the Grassmann leg ordering and the Array blocking.  A
single flat list of generators gives one binary dimension per generator.  A list
of groups gives one super-index per group.  A group with `n` generators becomes
an Array dimension of length `2^n`.

By default `index_order=:parity`, so each grouped index is ordered by even-parity
bitstrings followed by odd-parity bitstrings.  Within each parity block the
bitstrings are lexicographic in the generator order.  For two bits this is
`(0,0)=>1`, `(1,1)=>2`, `(0,1)=>3`, `(1,0)=>4`.

Set `index_order=:lexicographic` or `index_order=:little_endian` only for
low-level debugging or compatibility with older binary-order arrays.

`substitutions` is forwarded to `Symbolics.substitute` coefficient-wise before
entries are read.  Use, for example, `Dict(model.params.M1 => 1.0)` or the
`parameter_substitutions` helper.
"""
function coefficient_array(
    tensor::GrassmannExpr,
    leg_groups;
    substitutions=nothing,
    simplify::Bool=false,
    element_type=Any,
    index_order::Symbol=:parity,
)
    groups = _normalize_leg_groups(leg_groups)
    _validate_leg_groups(groups)
    prepared = _prepare_tensor_expression(tensor, substitutions)
    array = Array{element_type}(undef, coefficient_array_shape(groups))
    return _fill_coefficient_array!(
        array,
        bits -> get_coeff(prepared, occupied_legs_from_bits(groups, bits)),
        groups;
        simplify=simplify,
        index_order=index_order,
    )
end

"""
    coefficient_array(coefficient_at_bits::Function, leg_groups; substitutions=nothing,
                      simplify=false, element_type=Any, index_order=:parity)

Materialize an Array from a coefficient oracle.  The oracle is called as
`coefficient_at_bits(bits)`, where `bits` is the flattened occupation-bit tuple
in the order implied by `leg_groups`.  The oracle may return either a scalar or a
scalar `GrassmannExpr`.
"""
function coefficient_array(
    coefficient_at_bits::Function,
    leg_groups;
    substitutions=nothing,
    simplify::Bool=false,
    element_type=Any,
    index_order::Symbol=:parity,
)
    widths = _leg_group_widths(leg_groups)
    array = Array{element_type}(undef, Tuple(2^width for width in widths))
    return _fill_coefficient_array!(
        array,
        coefficient_at_bits,
        widths;
        substitutions=substitutions,
        simplify=simplify,
        index_order=index_order,
    )
end

const coefficient_tensor = coefficient_array

function parameter_substitutions(parameters::NamedTuple, values::NamedTuple)
    rules = Dict{Any,Any}()
    for name in keys(values)
        haskey(parameters, name) || throw(ArgumentError("unknown parameter name: $(name)"))
        rules[getproperty(parameters, name)] = getproperty(values, name)
    end
    return rules
end

parameter_substitutions(parameters::NamedTuple, values::AbstractDict) =
    Dict(getproperty(parameters, Symbol(name)) => value for (name, value) in values)