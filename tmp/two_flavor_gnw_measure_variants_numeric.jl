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
function vars_with_measure(kind)
    vars=_tfgnw_generators(); psi=vars.psi; psibar=vars.psibar
    mo = if kind == :current
        vars.measure_order
    elseif kind == :old_interleaved
        GrassmannGenerator[psi[1,1],psibar[1,1],psi[1,2],psibar[1,2],psi[2,1],psibar[2,1],psi[2,2],psibar[2,2]]
    elseif kind == :bar_then_psi_each_flavor
        GrassmannGenerator[psibar[1,1],psibar[1,2],psi[1,1],psi[1,2],psibar[2,1],psibar[2,2],psi[2,1],psi[2,2]]
    elseif kind == :all_psi_all_psibar
        GrassmannGenerator[psi[1,1],psi[1,2],psi[2,1],psi[2,2],psibar[1,1],psibar[1,2],psibar[2,1],psibar[2,2]]
    else
        error(kind)
    end
    merge(vars, (;measure_order=mo))
end
function numeric_pure(kind; mu1=0.0, mu2=0.0)
    vars=vars_with_measure(kind)
    params=(;M1=0.37,M2=-0.21,gs2=0.44,gp2=0.13,H=0.17,q=1/sqrt(2),z1=exp(mu1/2),z2=exp(mu2/2))
    base=_tfgnw_base_integrand(vars, params)
    pure=_tfgnw_integrate(vars, base; simplify=false)
    (;vars,pure)
end
function legs_direct(vars, idx)
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
function make_array(model)
    out=zeros(ComplexF64,16,16,16,16)
    for I in CartesianIndices(out)
        out[I]=ComplexF64(get_coeff(model.pure, legs_direct(model.vars,Tuple(I))))
    end
    out
end
function summarize(kind, ours, ext)
    maxdiff=0.0; maxidx=CartesianIndex(1,1,1,1); ndiff=0; both=0; onlyo=0; onlye=0; plus=0; minus=0; other=0
    for I in CartesianIndices(ours)
        o=ours[I]; e=ext[I]; d=abs(o-e)
        if d>maxdiff maxdiff=d; maxidx=I end
        d>1e-10 && (ndiff+=1)
        oz=abs(o)<1e-10; ez=abs(e)<1e-10
        if !oz && !ez
            both+=1; rr=o/e
            if abs(rr-1)<1e-8 plus+=1 elseif abs(rr+1)<1e-8 minus+=1 else other+=1 end
        elseif !oz && ez onlyo+=1 elseif oz && !ez onlye+=1 end
    end
    println(kind, " maxdiff=",maxdiff," maxidx=",Tuple(maxidx)," ours=",ours[maxidx]," ext=",ext[maxidx]," ndiff=",ndiff," both=",both," onlyo=",onlyo," onlye=",onlye," plus=",plus," minus=",minus," other=",other)
end
ext=Gross_Neveu_Wilson_model_Nf2(0.0,0.0,0.37,-0.21,0.44,0.13,0.17)
for kind in (:current,:old_interleaved,:bar_then_psi_each_flavor,:all_psi_all_psibar)
    println("KIND ",kind); flush(stdout)
    summarize(kind, make_array(numeric_pure(kind)), ext); flush(stdout)
end
