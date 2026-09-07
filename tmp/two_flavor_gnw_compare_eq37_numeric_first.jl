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
function numeric_model(; mu1=0.0, mu2=0.0, m1=0.37, m2=-0.21, gs2=0.44, gp2=0.13, H=0.17)
    vars=_tfgnw_generators()
    params=(;M1=m1,M2=m2,gs2=gs2,gp2=gp2,H=H,q=1/sqrt(2),z1=exp(mu1/2),z2=exp(mu2/2))
    base=_tfgnw_base_integrand(vars, params)
    pure=_tfgnw_integrate(vars, base; simplify=false)
    n11=_tfgnw_density(vars,1,1); n12=_tfgnw_density(vars,1,2); n21=_tfgnw_density(vars,2,1); n22=_tfgnw_density(vars,2,2)
    s11=_tfgnw_integrate(vars, n11*base; simplify=false)
    s12=_tfgnw_integrate(vars, n12*base; simplify=false)
    s21=_tfgnw_integrate(vars, n21*base; simplify=false)
    s22=_tfgnw_integrate(vars, n22*base; simplify=false)
    chiral=s11+s12+s21+s22
    singlet=im*(s11-s12+s21-s22)
    triplet=im*(s11-s12-s21+s22)
    (;vars,pure,chiral,singlet,triplet)
end
function legs_eq37(vars, idx)
    b1=F4[idx[1]]; b2=F4[idx[2]]; b3=F4[idx[3]]; b4=F4[idx[4]]
    i1,j1,k1,l1=b1; i2,j2,k2,l2=b2; i1p,j1p,k1p,l1p=b3; i2p,j2p,k2p,l2p=b4
    occ=(i1,j1,k1,l1, i2,j2,k2,l2, l1p,k1p,j1p,i1p, l2p,k2p,j2p,i2p)
    GrassmannGenerator[g for (g,b) in zip(vars.legs, occ) if b==1]
end
function make_array(tensor, vars)
    out=zeros(ComplexF64,16,16,16,16)
    for I in CartesianIndices(out)
        out[I]=ComplexF64(get_coeff(tensor, legs_eq37(vars,Tuple(I))))
    end
    out
end
function summarize(name, ours, ext)
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
    println(name)
    println("  max_abs_diff=",maxdiff," maxidx=",Tuple(maxidx)," ours=",ours[maxidx]," ext=",ext[maxidx])
    println("  nonzero_diff_count=",ndiff," support_both=",both," only_ours=",onlyo," only_ext=",onlye," plus=",plus," minus=",minus," other=",other)
end
for (label,mu1,mu2) in (("mu_zero",0.0,0.0),("mu_nonzero",0.31,-0.23))
    println("PARAM ",label); flush(stdout)
    model=numeric_model(;mu1,mu2)
    println("terms pure/chiral/singlet/triplet=",length(model.pure.terms),"/",length(model.chiral.terms),"/",length(model.singlet.terms),"/",length(model.triplet.terms)); flush(stdout)
    summarize("pure_eq37_numeric_first", make_array(model.pure,model.vars), Gross_Neveu_Wilson_model_Nf2(mu1,mu2,0.37,-0.21,0.44,0.13,0.17))
    if label == "mu_zero"
        ch,ps,pt=Gross_Neveu_Wilson_model_Nf2_impurity(0.37,-0.21,0.44,0.13,0.17)
        summarize("chiral_eq37_numeric_first", make_array(model.chiral,model.vars), ch)
        summarize("pseudoscalar_singlet_eq37_numeric_first", make_array(model.singlet,model.vars), ps)
        summarize("pseudoscalar_triplet_eq37_numeric_first", make_array(model.triplet,model.vars), pt)
    end
    flush(stdout)
end
