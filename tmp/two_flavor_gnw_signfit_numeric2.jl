using Pkg
const PROJECT_ROOT = dirname(@__DIR__)
Pkg.activate(PROJECT_ROOT)
using GrassmannSymbolics
include(joinpath(PROJECT_ROOT, "examples", "two_flavor_Gross_Neveu_Wilson.jl"))
include(get(ENV, "GRASSMANNTN_TENSOR_JL", joinpath(PROJECT_ROOT, "external_tensor.jl")))
const F4 = [
    (0,0,0,0), (1,1,0,0), (1,0,1,0), (1,0,0,1),
    (0,1,1,0), (0,1,0,1), (0,0,1,1), (1,1,1,1),
    (1,0,0,0), (0,1,0,0), (0,0,1,0), (1,1,1,0),
    (0,0,0,1), (1,1,0,1), (1,0,1,1), (0,1,1,1),
]
function numeric_pure()
    vars=_tfgnw_generators()
    params=(;M1=0.37,M2=-0.21,gs2=0.44,gp2=0.13,H=0.17,q=1/sqrt(2),z1=1.0,z2=1.0)
    pure=_tfgnw_integrate(vars, _tfgnw_base_integrand(vars, params); simplify=false)
    (;vars,pure)
end
function legs_direct(vars, idx)
    b1=F4[idx[1]]; b2=F4[idx[2]]; b3=F4[idx[3]]; b4=F4[idx[4]]
    i1,j1,k1,l1=b1; i2,j2,k2,l2=b2; i1p,j1p,k1p,l1p=b3; i2p,j2p,k2p,l2p=b4
    ordered=GrassmannGenerator[vars.alpha[1],vars.beta[1],vars.alpha[2],vars.beta[2], vars.eta[1],vars.zeta[1],vars.eta[2],vars.zeta[2], vars.betabar[1],vars.alphabar[1],vars.betabar[2],vars.alphabar[2], vars.zetabar[1],vars.etabar[1],vars.zetabar[2],vars.etabar[2]]
    occ=(i1,j1,k1,l1, i2,j2,k2,l2, j1p,i1p,l1p,k1p, j2p,i2p,l2p,k2p)
    GrassmannGenerator[g for (g,b) in zip(ordered,occ) if b==1]
end
bits_tuple(idx)=Tuple(vcat(collect(F4[idx[1]]),collect(F4[idx[2]]),collect(F4[idx[3]]),collect(F4[idx[4]])))
function gf2_solve(mat,b)
    m=length(mat); n=length(mat[1]); M=[copy(row) for row in mat]; rhs=copy(b); piv=Int[]; row=1
    for col in 1:n
        p=0
        for r in row:m
            if M[r][col]==1 p=r; break end
        end
        p==0 && continue
        M[row],M[p]=M[p],M[row]; rhs[row],rhs[p]=rhs[p],rhs[row]
        for r in 1:m
            if r!=row && M[r][col]==1
                @inbounds for c in col:n M[r][c] ⊻= M[row][c] end
                rhs[r] ⊻= rhs[row]
            end
        end
        push!(piv,col); row+=1; row>m && break
    end
    for r in row:m
        if all(==(0),M[r]) && rhs[r]==1 return nothing,piv end
    end
    x=zeros(Int,n); for (r,c) in enumerate(piv) x[c]=rhs[r] end
    x,piv
end
function feat(bits)
    xs=collect(bits); vals=Int[1]; append!(vals,xs)
    for mat in 1:16, b in mat+1:16 push!(vals,xs[mat]&xs[b]) end
    vals
end
function names()
    ns=String["1"]; append!(ns,["b$i" for i in 1:16]); for mat in 1:16, b in mat+1:16 push!(ns,"b$mat*b$b") end; ns
end
model=numeric_pure(); ext=Gross_Neveu_Wilson_model_Nf2(0.0,0.0,0.37,-0.21,0.44,0.13,0.17)
mat=Vector{Vector{Int}}(); y=Int[]
for I in CartesianIndices(ext)
    idx=Tuple(I); o=ComplexF64(get_coeff(model.pure, legs_direct(model.vars,idx))); e=ext[I]
    if abs(o)>1e-10 && abs(e)>1e-10
        push!(mat, feat(bits_tuple(idx)))
        push!(y, abs(o/e + 1)<1e-8 ? 1 : 0)
    end
end
x,piv=gf2_solve(mat,y)
println("eqs=",length(mat)," features=",length(mat[1])," solvable=",x!==nothing," rank=",length(piv))
if x!==nothing
    ns=names(); used=findall(==(1),x)
    println("terms_count=",length(used))
    println(join(ns[used]," + "))
    mism=0
    for i in eachindex(mat)
        pred=0; for j in used pred ⊻= mat[i][j] end
        pred==y[i] || (mism+=1)
    end
    println("mismatches=",mism)
end
