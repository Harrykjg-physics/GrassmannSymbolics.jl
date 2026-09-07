using Pkg
const PROJECT_ROOT = dirname(@__DIR__)
Pkg.activate(PROJECT_ROOT)
using GrassmannSymbolics
using Symbolics

include(joinpath(PROJECT_ROOT, "examples", "two_flavor_Gross_Neveu_Wilson.jl"))
external_tensor_path = get(ENV, "GRASSMANNTN_TENSOR_JL", joinpath(PROJECT_ROOT, "external_tensor.jl"))
isfile(external_tensor_path) || error("external tensor.jl not found: ", external_tensor_path)
include(external_tensor_path)

const GRASSMANNTN_F4_BITS = [
    (0,0,0,0), (1,1,0,0), (1,0,1,0), (1,0,0,1),
    (0,1,1,0), (0,1,0,1), (0,0,1,1), (1,1,1,1),
    (1,0,0,0), (0,1,0,0), (0,0,1,0), (1,1,1,0),
    (0,0,0,1), (1,1,0,1), (1,0,1,1), (0,1,1,1),
]

function derive_two_flavor_gnw_pure_only()
    @variables M1 M2 gs2 gp2 H q z1 z2
    params = (; M1, M2, gs2, gp2, H, q, z1, z2)
    vars = _tfgnw_generators()
    pure = _tfgnw_integrate(vars, _tfgnw_base_integrand(vars, params))
    return (; pure, params, vars)
end

function rules(model; mu1=0.0, mu2=0.0, m1=0.37, m2=-0.21, gs2=0.44, gp2=0.13, H=0.17)
    p = model.params
    Dict(p.M1=>m1, p.M2=>m2, p.gs2=>gs2, p.gp2=>gp2, p.H=>H, p.q=>1/sqrt(2), p.z1=>exp(mu1/2), p.z2=>exp(mu2/2))
end

function numeric_value(x)
    y = try Symbolics.simplify(x; expand=true) catch; x end
    try
        return ComplexF64(y)
    catch
        return Core.eval(Main, Meta.parse(string(y)))
    end
end

function occupied_legs(vars, index_tuple)
    b1 = GRASSMANNTN_F4_BITS[index_tuple[1]]
    b2 = GRASSMANNTN_F4_BITS[index_tuple[2]]
    b3 = GRASSMANNTN_F4_BITS[index_tuple[3]]
    b4 = GRASSMANNTN_F4_BITS[index_tuple[4]]
    i1,j1,k1,l1 = b1
    i2,j2,k2,l2 = b2
    i1p,j1p,k1p,l1p = b3
    i2p,j2p,k2p,l2p = b4
    ordered = GrassmannGenerator[
        vars.alpha[1], vars.beta[1], vars.alpha[2], vars.beta[2],
        vars.eta[1], vars.zeta[1], vars.eta[2], vars.zeta[2],
        vars.betabar[1], vars.alphabar[1], vars.betabar[2], vars.alphabar[2],
        vars.zetabar[1], vars.etabar[1], vars.zetabar[2], vars.etabar[2],
    ]
    occ = (i1,j1,k1,l1, i2,j2,k2,l2, j1p,i1p,l1p,k1p, j2p,i2p,l2p,k2p)
    return GrassmannGenerator[g for (g, bit) in zip(ordered, occ) if bit == 1]
end

function ours_array(model, r)
    prepared = substitute_coefficients(model.pure, r)
    out = Array{ComplexF64}(undef, 16,16,16,16)
    for I in CartesianIndices(out)
        out[I] = ComplexF64(numeric_value(get_coeff(prepared, occupied_legs(model.vars, Tuple(I)))))
    end
    out
end

function classify_ratio(z; tol=1e-8)
    candidates = Dict("+1"=>1+0im, "-1"=>-1+0im, "+i"=>0+1im, "-i"=>0-1im)
    best = nothing
    bestd = Inf
    for (name, val) in candidates
        d = abs(z-val)
        if d < bestd
            best = name; bestd = d
        end
    end
    return bestd < tol ? best : "other"
end

function diagnose(label; mu1=0.0, mu2=0.0)
    println("DIAG ", label)
    model = derive_two_flavor_gnw_pure_only()
    println("pure_terms=", length(model.pure.terms))
    r = rules(model; mu1, mu2)
    ours = ours_array(model, r)
    ext = Gross_Neveu_Wilson_model_Nf2(mu1, mu2, 0.37, -0.21, 0.44, 0.13, 0.17)
    both = 0; only_ours = 0; only_ext = 0; both_zero = 0
    counts = Dict("+1"=>0, "-1"=>0, "+i"=>0, "-i"=>0, "other"=>0)
    examples = []
    maxdiff = 0.0; maxidx = CartesianIndex(1,1,1,1)
    for I in CartesianIndices(ours)
        o = ours[I]; e = ext[I]
        d = abs(o-e)
        if d > maxdiff
            maxdiff = d; maxidx = I
        end
        oz = abs(o) < 1e-10; ez = abs(e) < 1e-10
        if oz && ez
            both_zero += 1
        elseif !oz && ez
            only_ours += 1
            length(examples) < 12 && push!(examples, (Tuple(I), "only_ours", o, e, missing))
        elseif oz && !ez
            only_ext += 1
            length(examples) < 12 && push!(examples, (Tuple(I), "only_ext", o, e, missing))
        else
            both += 1
            ratio = o/e
            cls = classify_ratio(ratio)
            counts[cls] = counts[cls] + 1
            if cls != "+1" && length(examples) < 12
                push!(examples, (Tuple(I), cls, o, e, ratio))
            end
        end
    end
    println("maxdiff=", maxdiff, " maxidx=", Tuple(maxidx), " ours=", ours[maxidx], " ext=", ext[maxidx])
    println("support both=", both, " only_ours=", only_ours, " only_ext=", only_ext, " both_zero=", both_zero)
    println("ratio_counts=", counts)
    for ex in examples
        idx, cls, o, e, ratio = ex
        bits = map(i -> GRASSMANNTN_F4_BITS[i], idx)
        println("example idx=", idx, " bits=", bits, " class=", cls, " ours=", o, " ext=", e, " ratio=", ratio)
    end
end

diagnose("mu_zero"; mu1=0.0, mu2=0.0)
diagnose("mu_nonzero"; mu1=0.31, mu2=-0.23)