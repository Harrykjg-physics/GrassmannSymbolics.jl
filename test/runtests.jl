using Test
using GrassmannSymbolics
using Symbolics

# ─────────────────────────────────────────────────────────────────────────────
# Helper: extract the coefficient of a monomial from a GrassmannExpr,
# falling back to 0 if the monomial is absent.
coeff(e::GrassmannExpr, gens) = get_coeff(e, collect(gens))
coeff(e::GrassmannExpr) = scalar_part(e)   # scalar part

# ─────────────────────────────────────────────────────────────────────────────
@testset "GrassmannSymbolics" begin

    # ── 1. Generator construction & ordering ─────────────────────────────────
    @testset "Generators" begin
        η = grassmann(:η)
        θ = grassmann(:θ)
        ηbar = dual(η)
        η1, η2, η3 = grassmann(:η, 1:3)

        @test η.name == :η
        @test η.index === nothing
        @test η.bar == false

        @test ηbar.bar == true
        @test ηbar.name == :η
        @test dual(ηbar) == η     # double dual = identity

        @test η1.index == 1
        @test η2.index == 2
        @test η3.index == 3

        # canonical ordering: η < η̄ < η₁ < η₂ < θ < θ̄
        @test η < ηbar
        @test η < η1
        @test η1 < η2
        @test η < θ
        @test ηbar < η1
        @test η2 < θ

        # @grassmann macro
        @grassmann α β γ
        @test α == GrassmannGenerator(:α, nothing, false)
        @test β.name == :β
        @test γ.name == :γ
    end

    # ── 2. GrassmannExpr construction ────────────────────────────────────────
    @testset "Construction" begin
        @grassmann η θ

        e1 = GrassmannExpr(η)
        @test !iszero(e1)
        @test coeff(e1, [η]) == 1

        e_scalar = GrassmannExpr(5)
        @test scalar_part(e_scalar) == 5
        @test !is_grassmann(e_scalar)

        e_zero = GrassmannExpr(0)
        @test iszero(e_zero)
    end

    # ── 3. Addition & subtraction ─────────────────────────────────────────────
    @testset "Addition / Subtraction" begin
        @grassmann η θ

        s = η + θ
        @test coeff(s, [η]) == 1
        @test coeff(s, [θ]) == 1
        @test coeff(s) == 0

        s2 = 2η + 3θ
        @test coeff(s2, [η]) == 2
        @test coeff(s2, [θ]) == 3

        @test iszero(η - η)
        @test coeff(η + 1, [η]) == 1
        @test scalar_part(η + 1) == 1

        @test coeff(-η, [η]) == -1
        @test coeff(η - θ, [θ]) == -1
    end

    # ── 4. Scalar multiplication / division ──────────────────────────────────
    @testset "Scalar multiply / divide" begin
        @grassmann η θ

        @test coeff(3η, [η]) == 3
        @test coeff(η * 3, [η]) == 3
        @test coeff(η / 2, [η]) == 1 // 2 || coeff(η / 2, [η]) ≈ 0.5
        @test iszero(0 * η)
    end

    # ── 5. Anticommutative multiplication ────────────────────────────────────
    @testset "Exterior product (anticommutativity)" begin
        @grassmann η θ ψ

        eηθ = η * θ
        eθη = θ * η

        @test coeff(eηθ, [η, θ]) == 1
        @test coeff(eθη, [η, θ]) == -1   # θη = −ηθ

        # The product ηθ as a grade-2 element commutes with any other element
        # (even parity: (ηθ)·ψ = ψ·(ηθ))
        lhs = (η * θ) * ψ
        rhs = ψ * (η * θ)
        @test coeff(lhs, [η, θ, ψ]) == coeff(rhs, [η, θ, ψ])

        # Triple product: ηθψ should appear with sign accounting for permutation
        @test coeff(η * θ * ψ, [η, θ, ψ]) == 1
        @test coeff(θ * η * ψ, [η, θ, ψ]) == -1  # one swap
        @test coeff(θ * ψ * η, [η, θ, ψ]) == 1   # two swaps

        # More products
        @test coeff(η * θ + θ * η, [η, θ]) == 0   # η θ + θ η = 0

        # Nilpotency
        @test iszero(η * η)
        @test iszero(η^2)
        @test iszero((η + θ)^3)   # at most 2 generators involved ⟹ grade ≤ 2
    end

    # ── 6. Power ──────────────────────────────────────────────────────────────
    @testset "Power" begin
        @grassmann η θ

        @test coeff(η^0) == 1
        @test iszero(grassmann_part(η^0))

        @test coeff(η^1, [η]) == 1
        @test iszero(η^2)
        @test iszero(η^3)

        eηθ = η * θ
        @test coeff(eηθ^0) == 1
        @test coeff(eηθ^1, [η, θ]) == 1
        @test iszero(eηθ^2)   # (ηθ)² = ηθηθ = -ηηθθ = 0
    end

    # ── 7. exp ────────────────────────────────────────────────────────────────
    @testset "exp" begin
        @grassmann η θ

        # exp(η) = 1 + η
        e1 = exp(η)
        @test scalar_part(e1) == 1
        @test coeff(e1, [η]) == 1
        @test length(e1.terms) == 2

        # exp(ηθ) = 1 + ηθ
        e2 = exp(η * θ)
        @test scalar_part(e2) == 1
        @test coeff(e2, [η, θ]) == 1

        # exp(c + η) = exp(c) * (1 + η)  for numeric c
        c = 2.0
        e3 = exp(GrassmannExpr(c) + GrassmannExpr(η))
        @test scalar_part(e3) ≈ exp(c)
        @test coeff(e3, [η]) ≈ exp(c)

        # exp(0) = 1
        @test scalar_part(exp(GrassmannExpr(0))) == 1

        # Nilpotency: exp truncates automatically
        # Three generators: exp(ηθψ) = 1 + ηθψ  (grade > 3 = 0)
        @grassmann ψ
        e4 = exp(η * θ * ψ)
        @test scalar_part(e4) == 1
        @test coeff(e4, [η, θ, ψ]) == 1
    end

    # ── 8. Predicates ─────────────────────────────────────────────────────────
    @testset "Predicates" begin
        @grassmann η θ

        @test is_even(GrassmannExpr(1))         # scalar is even
        @test is_even(η * θ)                    # grade 2 is even
        @test !is_even(η)                       # grade 1 is odd
        @test is_odd(η)
        @test !is_odd(η * θ)
        @test is_grassmann(η)
        @test !is_grassmann(GrassmannExpr(5))
        @test is_grassmann(η + 1)               # mixed: has Grassmann part
    end

    # ── 9. get_coeff & generators ─────────────────────────────────────────────
    @testset "get_coeff / generators()" begin
        @grassmann η θ ψ

        e = 3(η * θ) + 2η + 5
        @test coeff(e, [η, θ]) == 3
        @test coeff(e, [θ, η]) == -3    # anticommutation sign
        @test coeff(e, [η]) == 2
        @test scalar_part(e) == 5
        @test coeff(e, [ψ]) == 0

        gs = generators(e)
        @test η ∈ gs
        @test θ ∈ gs
        @test ψ ∉ gs
    end

    # ── 10. Berezin integral ──────────────────────────────────────────────────
    @testset "Berezin integral" begin
        @grassmann η θ ψ

        # ∫dη η = 1
        @test scalar_part(d(η) * GrassmannExpr(η)) == 1

        # ∫dη 1 = 0
        @test iszero(d(η) * GrassmannExpr(1))

        # ∫dη θ = 0   (θ ≠ η)
        @test iszero(d(η) * GrassmannExpr(θ))

        # ∫dη (ηθ) = θ
        r1 = d(η) * (η * θ)
        @test coeff(r1, [θ]) == 1
        @test iszero(grassmann_part(r1) - GrassmannExpr(θ))

        # ∫dθ (ηθ) = η  (sign: θ is at position 2, sign = (−1)^1 = −1,
        # but canonical monomial is sorted as ηθ so θ is at index 2)
        r2 = d(θ) * (η * θ)
        @test coeff(r2, [η]) == -1  # -η

        # ∫dη dθ (ηθ) — dθ applied first:
        #   ∫dθ(ηθ) = -η
        #   ∫dη(-η) = -1
        r3 = d(η) * d(θ) * (η * θ)
        @test scalar_part(r3) == -1

        # ∫dθ dη (ηθ) — dη applied first:
        #   ∫dη(ηθ) = θ
        #   ∫dθ(θ)  = 1
        r4 = d(θ) * d(η) * (η * θ)
        @test scalar_part(r4) == 1

        # Combining measures with *
        meas = d(θ) * d(η)   # BerezinMeasure
        @test meas isa BerezinMeasure
        r5 = meas * (η * θ)
        @test scalar_part(r5) == 1

        # Repeated measure → 0
        @test iszero(d(η) * d(η) * GrassmannExpr(η))

        # ∫dη dθ dψ (ηθψ) — should give ±1
        r6 = d(η) * d(θ) * d(ψ) * (η * θ * ψ)
        @test abs(scalar_part(r6)) == 1

        # Integration of a sum
        e = η * θ + η * ψ
        r7 = d(η) * e
        @test coeff(r7, [θ]) == 1
        @test coeff(r7, [ψ]) == 1
    end

    # ── 11. Display (smoke-test that show doesn't throw) ─────────────────────
    @testset "Display" begin
        @grassmann η θ
        ηbar = dual(η)
        η1, η2 = grassmann(:η, 1:2)

        @test sprint(show, η) == "η"
        @test sprint(show, ηbar) isa String   # contains overbar
        @test sprint(show, η1) == "η₁"
        @test sprint(show, η2) == "η₂"

        e = 2η + 3(η * θ)
        s = sprint(show, e)
        @test occursin("η", s)
        @test occursin("θ", s)

        @test sprint(show, GrassmannExpr(0)) == "0"
        @test sprint(show, d(η)) isa String
    end

    # ── 12. Dual generators ───────────────────────────────────────────────────
    @testset "Dual generators" begin
        @grassmann η θ
        ηbar = dual(η)
        θbar = dual(θ)

        # η and η̄ anticommute
        p = η * ηbar
        q = ηbar * η
        @test coeff(p, [η, ηbar]) == 1
        @test coeff(q, [η, ηbar]) == -1

        # (η + η̄)² = ηη̄ + η̄η = 0  (because ηη̄ = -η̄η)
        @test iszero((η + ηbar)^2)

        # Cross terms: η * θ̄  etc.
        p2 = η * θbar + θbar * η
        @test iszero(p2)   # anticommutation again
    end

    # ── 13. Mixed scalar + Grassmann ─────────────────────────────────────────
    @testset "Mixed expressions" begin
        @grassmann η θ

        e = 1 + η + η * θ
        @test scalar_part(e) == 1
        @test coeff(e, [η]) == 1
        @test coeff(e, [η, θ]) == 1

        # (1 + η)(1 + θ) = 1 + θ + η + ηθ
        prod = (1 + GrassmannExpr(η)) * (1 + GrassmannExpr(θ))
        @test scalar_part(prod) == 1
        @test coeff(prod, [η]) == 1
        @test coeff(prod, [θ]) == 1
        @test coeff(prod, [η, θ]) == 1

        # (1 + η)² = 1 + 2η + η² = 1 + 2η  (η² = 0)
        sq = (1 + GrassmannExpr(η))^2
        @test scalar_part(sq) == 1
        @test coeff(sq, [η]) == 2
    end

    # ── 14. Symbolic utilities and local tensor kernel ────────────────────────
    @testset "Coefficient utilities / local tensor" begin
        @grassmann ψ η
        ψbar = dual(ψ)
        ηbar = dual(η)

        @test -η == (-1) * η
        @test grassmann(:ξ, 1:2; bar=true) == dual.(grassmann(:ξ, 1:2))

        expression = 2η + 3(η * ηbar)
        doubled = map_coefficients(c -> 2c, expression)
        @test get_coeff(doubled, [η]) == 4
        @test get_coeff(doubled, [η, ηbar]) == 6
        @test symbolically_equal(doubled, 2 * expression)

        @variables a
        substituted = substitute_coefficients(a * η, Dict(a => 7))
        @test get_coeff(substituted, [η]) == 7

        @grassmann θ
        renamed = substitute_generators(η * θ, Dict(η => θ, θ => η))
        @test get_coeff(renamed, [η, θ]) == -1
        @test iszero(map_generators(_ -> η, η * θ))

        # ∫dψ dψbar exp(-3ψbarψ)
        #   exp(-ψbarη + 1/2 ηbarψ) = -3 - 1/2 ηηbar.
        tensor = local_grassmann_tensor(
            3 * ψbar * ψ,
            [-ψbar * η + (1 // 2) * ηbar * ψ],
            [ψ, ψbar],
        )
        @test scalar_part(tensor) == -3
        @test get_coeff(tensor, [η, ηbar]) == -(1 // 2)
        @test is_even(tensor)

        @grassmann chi eta
        chibar = dual(chi)
        etabar = dual(eta)
        ascii_tensor = local_grassmann_tensor(
            3 * chibar * chi,
            [-chibar * eta + (1 // 2) * etabar * chi],
            [chi, chibar],
        )
        channels = split_hopping_channel(
            chibar,
            eta,
            chi;
            left_scale=-1,
            right_scale=1 // 2,
            tag=:test_hopping,
        )
        spec = NearestNeighborTensorSpec(
            3 * chibar * chi,
            [chi, chibar],
            channels;
            leg_order=[eta, etabar],
            metadata=Dict(:model => :one_component_test),
        )
        compiled = compile_nearest_neighbor_tensor(spec)
        @test symbolically_equal(compiled.tensor, ascii_tensor)
        @test compiled.legs == [eta, etabar]
        @test compiled.metadata[:model] == :one_component_test
        @test compiled.channels[1].tag == :test_hopping

        action = NearestNeighborAction(
            [chi, chibar];
            terms=[
                OnsiteTerm(3 * chibar * chi; tag=:mass),
                FactorizedHoppingTerm(
                    chibar,
                    eta,
                    chi;
                    left_scale=-1,
                    right_scale=1 // 2,
                    tag=:test_hopping,
                ),
            ],
            leg_order=[eta, etabar],
            metadata=Dict(:model => :one_component_action_test),
        )
        action_spec = nearest_neighbor_tensor_spec(action)
        action_compiled = compile_nearest_neighbor_tensor(action)
        @test symbolically_equal(action_compiled.tensor, ascii_tensor)
        @test action_compiled.legs == [eta, etabar]
        @test length(action_spec.channels) == 2
        @test action_compiled.metadata[:model] == :one_component_action_test
        @test action_compiled.metadata[:term_tags] == [:mass, :test_hopping]

        closed = contract_grassmann((1 + eta,), ((etabar, eta),))
        @test scalar_part(closed) == 1
        @test isempty(generators(closed))
    end

    @testset "Coefficient arrays" begin
        @grassmann a b c d
        expression = 5 + 2a + 3b * c - 7a * d
        groups = ((a, b), (c, d))

        @test occupation_bits_from_index(1, 2) == (0, 0)
        @test occupation_bits_from_index(2, 2) == (1, 1)
        @test occupation_bits_from_index(3, 2) == (0, 1)
        @test occupation_bits_from_index(4, 2) == (1, 0)
        @test occupation_index_from_bits((1, 1)) == 2
        @test occupation_index_from_bits((0, 1)) == 3
        @test occupation_index_from_bits((1, 0)) == 4

        array = coefficient_array(expression, groups)
        @test size(array) == (4, 4)
        @test array[1, 1] == 5
        @test array[4, 1] == 2
        @test array[3, 4] == 3
        @test array[4, 3] == -7

        binary_array = coefficient_array(expression, groups; index_order=:little_endian)
        @test binary_array[2, 1] == 2
        @test binary_array[3, 2] == 3
        @test binary_array[2, 3] == -7

        @variables x
        substituted = coefficient_array(x * a * b, ((a, b),); substitutions=Dict(x => 11))
        @test size(substituted) == (4,)
        @test substituted[2] == 11

        oracle_array = coefficient_array(bits -> sum(bits), groups)
        @test oracle_array[1, 1] == 0
        @test oracle_array[2, 2] == 4

        width_array = coefficient_array(bits -> sum(bits), (2, 2))
        @test size(width_array) == (4, 4)
        @test width_array[2, 2] == 4
    end
end # @testset "GrassmannSymbolics"

@testset "3D staggered fermion SMG example" begin
    example_path = joinpath(@__DIR__, "..", "examples", "3d_staggered_fermion_SMG.jl")
    @test isfile(example_path)

    if isfile(example_path)
        include(example_path)
        model = derive_3d_staggered_fermion_smg_tensor()
        p = model.parameters

        empty_bits = smg_bits(ntuple(_ -> 0, 36))
        @test symbolically_equal(
            GrassmannExpr(smg_local_integral(empty_bits, model)),
            GrassmannExpr(p.UI),
        )
        @test is_even(smg_tensor_entry(empty_bits, model))

        onsite_u_dimer = smg_bits((0, 0, 1, 0, 0, 0, ntuple(_ -> 0, 30)...))
        @test smg_local_integral(onsite_u_dimer, model) == 0

        single_u_hop = smg_bits((1, 0, 0, 0, 0, 0, ntuple(_ -> 0, 30)...))
        @test smg_local_integral(single_u_hop, model) == 0

        complementary_dimers = smg_bits((0, 0, 1, 0, 0, 1, ntuple(_ -> 0, 30)...))
        @test symbolically_equal(
            GrassmannExpr(smg_local_integral(complementary_dimers, model)),
            GrassmannExpr(p.sqrt_UB^2),
        )

        repeated_u_dimer = smg_bits((0, 0, 1, 0, 0, 0,
                                     0, 0, 0, 0, 0, 0,
                                     0, 0, 0, 0, 0, 0,
                                     0, 0, 1, 0, 0, 0,
                                     0, 0, 0, 0, 0, 0,
                                     0, 0, 0, 0, 0, 0))
        @test symbolically_equal(
            GrassmannExpr(smg_local_integral(repeated_u_dimer, model)),
            GrassmannExpr(0),
        )
    end
end

@testset "3D staggered fermion SMG no-UB example" begin
    example_path = joinpath(@__DIR__, "..", "examples", "3d_staggered_fermion_SMG_no_UB.jl")
    @test isfile(example_path)

    if isfile(example_path)
        include(example_path)
        model = derive_3d_staggered_fermion_smg_no_ub_tensor()
        p = model.parameters

        empty_bits = smg_no_ub_bits(ntuple(_ -> 0, 24))
        @test symbolically_equal(
            GrassmannExpr(smg_no_ub_local_integral(empty_bits, model)),
            GrassmannExpr(p.UI),
        )
        @test is_even(smg_no_ub_tensor_entry(empty_bits, model))
        @test coefficient_array_shape(model.leg_groups) == (16, 16, 16, 16, 16, 16)

        single_u_hop = smg_no_ub_bits((1, 0, 0, 0, ntuple(_ -> 0, 20)...))
        @test smg_no_ub_local_integral(single_u_hop, model) == 0

        saturated_hops = smg_no_ub_bits((1, 1, 1, 1, ntuple(_ -> 0, 20)...))
        @test symbolically_equal(
            GrassmannExpr(smg_no_ub_local_integral(saturated_hops, model)),
            GrassmannExpr(p.p_t^4),
        )
    end
end
