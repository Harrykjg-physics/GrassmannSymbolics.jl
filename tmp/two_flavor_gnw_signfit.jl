using Pkg
const PROJECT_ROOT = dirname(@__DIR__)
Pkg.activate(PROJECT_ROOT)
using GrassmannSymbolics
using Symbolics
include(joinpath(PROJECT_ROOT, "examples", "two_flavor_Gross_Neveu_Wilson.jl"))
external_tensor_path = get(ENV, "GRASSMANNTN_TENSOR_JL", joinpath(PROJECT_ROOT, "external_tensor.jl"))
include(external_tensor_path)
const F4 = [
    (0,0,0,0), (1,1,0,0), (1,0,1,0), (1,0,0,1),
    (0,1,1,0), (0,1,0,1), (0,0,1,1), (1,1,1,1),
    (1,0,0,0), (0,1,0,0), (0,0,1,0), (1,1,1,0),
    (0,0,0,1), (1,1,0,1), (1,0,1,1), (0,1,1,1),
]
function derive()
    @variables M1 M2 gs2 gp2 H q z1 z2
    params=(;M1,M2,gs2,gp2,H,q,z1,z2)
    vars=_tfgnw_generators()
    pure=_tfgnw_integrate(vars,_tfgnw_base_integrand(vars,params))
    (;pure,params,vars)
end
function rules(model)
    p=model.params
    Dict(p.M1=>0.37,p.M2=>-0.21,p.gs2=>0.44,p.gp2=>0.13,p.H=>0.17,p.q=>1/sqrt(2),p.z1=>1.0,p.z2=>1.0)
end
function numeric_value(x)
    y=try Symbolics.simplify(x; expand=true) catch; x end
    try ComplexF64(y) catch; Core.eval(Main, Meta.parse(string(y))) end
end
function legs(vars, idx)
    b1=F4[idx[1]]; b2=F4[idx[2]]; b3=F4[idx[3]]; b4=F4[idx[4]]
    i1,j1,k1,l1=b1; i2,j2,k2,l2=b2; i1p,j1p,k1p,l1p=b3; i2p,j2p,k2p,l2p=b4
    ordered=GrassmannGenerator[
        vars.alpha[1],vars.beta[1],vars.alpha[2],vars.beta[2],
        vars.eta[1],vars.zeta[1],vars.eta[2],vars.zeta[2],
        vars.betabar[1],vars.alphabar[1],vars.betabar[2],vars.alphabar[2],
        vars.zetabar[1],vars.etabar[1],vars.zetabar[2],vars.etabar[2],
    ]
    occ=(i1,j1,k1,l1, i2,j2,k2,l2, j1p,i1p,l1p,k1p, j2p,i2p,l2p,k2p)
    GrassmannGenerator[g for (g,b) in zip(ordered,occ) if b==1]
end
function bits_tuple(idx)
    b1=F4[idx[1]]; b2=F4[idx[2]]; b3=F4[idx[3]]; b4=F4[idx[4]]
    return Tuple(vcat(collect(b1), collect(b2), collect(b3), collect(b4)))
end
function gf2_solve(A,b)
    m=length(A); n=length(A[1]); M=[Vector{Int}(A[i]) for i in 1:m]; rhs=Vector{Int}(b)
    pivots=Int[]; row=1
    for col in 1:n
        p=findfirst(r -> M[r][col]==1, row:m)
        p===nothing && continue
        p = p + row - 1
        M[row],M[p]=M[p],M[row]; rhs[row],rhs[p]=rhs[p],rhs[row]
        for r in 1:m
            if r != row && M[r][col]==1
                M[r] = xor.(M[r], M[row]); rhs[r] ⊻= rhs[row]
            end
        end
        push!(pivots,col); row += 1
        row>m && break
    end
    for r in row:m
        if all(x->x==0,M[r]) && rhs[r]==1
            return nothing, pivots
        end
    end
    x=zeros(Int,n)
    for (r,col) in enumerate(pivots)
        x[col]=rhs[r]
    end
    return x,pivots
end
function feature(bits)
    xs=collect(bits)
    vals=Int[1]
    append!(vals,xs)
    for a in 1:16, b in a+1:16
        push!(vals, xs[a]*xs[b])
    end
    vals
end
function feature_names()
    names=String["1"]
    append!(names,["b$i" for i in 1:16])
    for a in 1:16, b in a+1:16
        push!(names,"b$a*b$b")
    end
    names
end
model=derive(); prep=substitute_coefficients(model.pure, rules(model)); ext=Gross_Neveu_Wilson_model_Nf2(0.0,0.0,0.37,-0.21,0.44,0.13,0.17)
A=Vector{Vector{Int}}(); b=Int[]; samples=[]
for idx in CartesianIndices(ext)
    t=Tuple(idx); o=ComplexF64(numeric_value(get_coeff(prep, legs(model.vars,t)))); e=ext[idx]
    if abs(o)>1e-10 && abs(e)>1e-10
        push!(A, feature(bits_tuple(t)))
        push!(b, abs(o/e + 1) < 1e-8 ? 1 : 0)
        length(samples)<10 && abs(o/e + 1)<1e-8 && push!(samples,(t,bits_tuple(t)))
    end
end
x,piv=gf2_solve(A,b)
println("equations=",length(A)," features=",length(A[1])," solvable=",x!==nothing," rank=",length(piv))
if x!==nothing
    names=feature_names(); used=findall(==(1),x)
    println("terms_count=",length(used))
    println(join(names[used], " + "))
    mism=0
    for i in eachindex(A)
        pred=foldl(⊻, (A[i][j] for j in used); init=0)
        pred==b[i] || (mism+=1)
    end
    println("fit_mismatches=",mism)
end
println("samples_minus=")
for s in samples println(s) end
