using Pkg
const PROJECT_ROOT = dirname(@__DIR__)
Pkg.activate(PROJECT_ROOT)
using GrassmannSymbolics
using Symbolics
include(joinpath(PROJECT_ROOT, "examples", "two_flavor_Gross_Neveu_Wilson.jl"))
include(get(ENV, "GRASSMANNTN_TENSOR_JL", joinpath(PROJECT_ROOT, "external_tensor.jl")))
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
    chiral=simplify_coefficients(_tfgnw_integrate(vars,_tfgnw_base_integrand(vars,params)*_tfgnw_chiral_operator(vars)))
    singlet=simplify_coefficients(_tfgnw_integrate(vars,_tfgnw_base_integrand(vars,params)*_tfgnw_pseudoscalar_singlet(vars)))
    triplet=simplify_coefficients(_tfgnw_integrate(vars,_tfgnw_base_integrand(vars,params)*_tfgnw_pseudoscalar_triplet(vars)))
    (;pure,chiral,singlet,triplet,params,vars)
end
function rules(model; mu1=0.0, mu2=0.0, m1=0.37, m2=-0.21, gs2=0.44, gp2=0.13, H=0.17)
    p=model.params
    Dict(p.M1=>m1,p.M2=>m2,p.gs2=>gs2,p.gp2=>gp2,p.H=>H,p.q=>1/sqrt(2),p.z1=>exp(mu1/2),p.z2=>exp(mu2/2))
end
function numeric_value(x)
    y=try Symbolics.simplify(x; expand=true) catch; x end
    try ComplexF64(y) catch; Core.eval(Main, Meta.parse(string(y))) end
end
function legs_eq37(vars, idx)
    b1=F4[idx[1]]; b2=F4[idx[2]]; b3=F4[idx[3]]; b4=F4[idx[4]]
    i1,j1,k1,l1=b1; i2,j2,k2,l2=b2; i1p,j1p,k1p,l1p=b3; i2p,j2p,k2p,l2p=b4
    occ=(i1,j1,k1,l1, i2,j2,k2,l2, l1p,k1p,j1p,i1p, l2p,k2p,j2p,i2p)
    GrassmannGenerator[g for (g,b) in zip(vars.legs, occ) if b==1]
end
function make_array(tensor, model, r)
    prep=substitute_coefficients(tensor,r)
    out=zeros(ComplexF64,16,16,16,16)
    for I in CartesianIndices(out)
        out[I]=ComplexF64(numeric_value(get_coeff(prep, legs_eq37(model.vars,Tuple(I)))))
    end
    out
end
function summarize(name, ours, ext)
    maxdiff=0.0; maxidx=CartesianIndex(1,1,1,1); nonzero=0; both=0; onlyo=0; onlye=0; ratios=Dict("+1"=>0,"-1"=>0,"other"=>0)
    for I in CartesianIndices(ours)
        d=abs(ours[I]-ext[I]); if d>maxdiff maxdiff=d; maxidx=I end
        d>1e-10 && (nonzero+=1)
        oz=abs(ours[I])<1e-10; ez=abs(ext[I])<1e-10
        if !oz && !ez
            both+=1; rr=ours[I]/ext[I]
            if abs(rr-1)<1e-8 ratios["+1"]+=1 elseif abs(rr+1)<1e-8 ratios["-1"]+=1 else ratios["other"]+=1 end
        elseif !oz && ez onlyo+=1 elseif oz && !ez onlye+=1 end
    end
    println(name)
    println("  max_abs_diff=",maxdiff," maxidx=",Tuple(maxidx)," ours=",ours[maxidx]," ext=",ext[maxidx])
    println("  nonzero_diff_count=",nonzero," support_both=",both," only_ours=",onlyo," only_ext=",onlye," ratios=",ratios)
end
model=derive(); println("terms pure/chiral/singlet/triplet=",length(model.pure.terms),"/",length(model.chiral.terms),"/",length(model.singlet.terms),"/",length(model.triplet.terms))
for (label,mu1,mu2) in (("mu_zero",0.0,0.0),("mu_nonzero",0.31,-0.23))
    println("PARAM ",label); r=rules(model;mu1,mu2)
    summarize("pure_eq37", make_array(model.pure,model,r), Gross_Neveu_Wilson_model_Nf2(mu1,mu2,0.37,-0.21,0.44,0.13,0.17))
end
println("IMPURITY mu_zero")
r=rules(model)
ch,ps,pt=Gross_Neveu_Wilson_model_Nf2_impurity(0.37,-0.21,0.44,0.13,0.17)
summarize("chiral_eq37", make_array(model.chiral,model,r), ch)
summarize("pseudoscalar_singlet_eq37", make_array(model.singlet,model,r), ps)
summarize("pseudoscalar_triplet_eq37", make_array(model.triplet,model,r), pt)
