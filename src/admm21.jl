using LinearAlgebra, JuMP,MosekTools, Gurobi
using LinearAlgebra
using JuMP, MosekTools, Gurobi
using MAT,DelimitedFiles, Combinatorics

import MutableArithmetics
const MA = MutableArithmetics
const SA = MA.SparseArrays

using JLD


include("types.jl")
include("util.jl")
include("probs.jl")
include("instances.jl")

function getnorm21(A::Matrix{T}) where {T}
    m,n = size(A);
    norm_21 = sum(norm(A[i,:],2) for i = (1:m));
    return norm_21
end

function min21norm(A)
    m,n = size(A);
    r = rank(A,1e-6);
    U,Σ,V = svd(A,full=true);
    D = Σ[1:r];
    Dinv = 1 ./D;
    Ur = U[:,1:r];
    Vr,Vs = V[:,1:r],V[:,r+1:n];
    G = Vr*diagm(Dinv)*Ur';
    model = Model(Mosek.Optimizer);
    set_time_limit_sec(model, 18000)
    set_silent(model);
    @variable(model, t[i in (1:n)]);
    @variable(model, Z[1:n-r,1:r]);
    Cmat = Vr*diagm(Dinv) + Vs*Z;
    @constraint(model, [i = (1:n)], [t[i];  Cmat[i,:]] in SecondOrderCone())
    @objective(model,Min,sum(t))
    optimize!(model);
    primal_status = termination_status(model);
    H = G +Vs*value.(Z)*Ur';
    z =  objective_value(model)
    return z#,H
end

function solve21msk(inst::GinvInst,ginvInit::GinvInit)
    n,r = inst.n,inst.r;
    model = Model(Mosek.Optimizer);
    set_time_limit_sec(model, 7200)
    set_silent(model);
    set_attribute(model, "BASIS_REL_TOL_S", 1e-5)
    set_attribute(model, "INTPNT_CO_TOL_DFEAS", 1e-5)
    set_attribute(model, "INTPNT_CO_TOL_INFEAS", 1e-5)
    set_attribute(model, "INTPNT_CO_TOL_MU_RED", 1e-5)
    set_attribute(model, "INTPNT_CO_TOL_PFEAS", 1e-5)
    set_attribute(model, "INTPNT_CO_TOL_REL_GAP", 1e-5)
    set_attribute(model, "INTPNT_TOL_DFEAS", 1e-5)
    set_silent(model);
    @variable(model, t[i in (1:n)]);
    @variable(model, Z[1:n-r,1:r]);
    Cmat = ginvInit.Vr*diagm(ginvInit.Dinv) + ginvInit.Vs*Z;
    @constraint(model, [i = (1:n)], [t[i];  Cmat[i,:]] in SecondOrderCone())
    @objective(model,Min,sum(t))
    time_to_solve = @elapsed optimize!(model);

    csol = SolutionOptimizer();
    primal_status = termination_status(model);
    println(string("Solution status Mosek: ", primal_status))
    if primal_status != MOI.TIME_LIMIT
        csol.H = ginvInit.G + ginvInit.Vs*value.(Z)*ginvInit.Ur';
        csol.z,csol.time =  objective_value(model),time_to_solve;
        csol.iter = 0;
    else
        csol.H = zeros(n,m);
        csol.z,csol.time = Inf, time_to_solve;
        csol.iter = 0;
    end
    return csol
end

function solve21grb(inst::GinvInst,ginvInit::GinvInit)
    n,r = inst.n,inst.r;
    model = Model(Gurobi.Optimizer);
    
    set_silent(model);
    set_time_limit_sec(model, 7200)
    set_attribute(model, "BarConvTol", 1e-5)
    set_attribute(model, "BarQCPConvTol", 1e-5)
    set_attribute(model, "FeasibilityTol", 1e-5)
    set_attribute(model, "OptimalityTol", 1e-5)
    set_silent(model);
    @variable(model, t[i in (1:n)]);
    @variable(model, Z[1:n-r,1:r]);
    Cmat = ginvInit.Vr*diagm(ginvInit.Dinv) + ginvInit.Vs*Z;
    @constraint(model, [i = (1:n)], [t[i];  Cmat[i,:]] in SecondOrderCone())
    @objective(model,Min,sum(t))
    time_to_solve = @elapsed optimize!(model);

    csol = SolutionOptimizer();
    primal_status = termination_status(model);
    println(string("Solution status Gurobi: ", primal_status))
    
    if primal_status != MOI.TIME_LIMIT
        csol.H = ginvInit.G + ginvInit.Vs*value.(Z)*ginvInit.Ur';
        csol.z,csol.time =  objective_value(model),time_to_solve;
        csol.iter = 0;
    else
        csol.H = zeros(n,m);
        csol.z,csol.time = Inf, time_to_solve;
        csol.iter = 0;
    end
    return csol
end

function update_E21n(V1DinvΛ,V2Z,mu)
    Y = V1DinvΛ + V2Z
    E = zeros(n,r)
    norm_rows_y = norm.(eachrow(Y))
    for (i,el) in enumerate(norm_rows_y)
        if (1/mu) < el
            E[i,:] = ((el - (1/mu))/el) * Y[i,:]
        end
    end
    return E
end

function update_Z21n(V1DinvΛ,V2,E)
    return V2'*(-V1DinvΛ + E)
end

function runADMM21n(V1Dinv,V2,U1,Λ)
    mu = 0.1
    #Λ = zeros(n,r)
    E = zeros(n,r)
    Z = zeros(n-r,r)
    V2Z = V2*Z;
    iter = 0
    while true
        V1DinvΛ = V1Dinv + (1/mu)*Λ
        E = update_E21n(V1DinvΛ,V2Z,mu)
        Z = update_Z21n(V1DinvΛ,V2,E)
        V2Z = V2*Z;
        res_infeas = V1Dinv + V2Z-E
        Λ += mu*res_infeas
        mu = 1.05*mu
        iter += 1
        if norm(res_infeas) < 1e-5
            break
        end
    end
    admmsol = SolutionOptimizer();
    admmsol.H = V1Dinv*U1' +V2*Z*U1';
    admmsol.iter = iter
    return admmsol
end

arr_res_all = [];
arr_resADMM,arr_resGRB,arr_resMSK = init_ginv_res(),init_ginv_res(),init_ginv_res();
global m,n,r  = 100,200,70;
M = [40,80,120,160,200,240,280,320,1000,2000,3000,4000,5000,6000,7000,8000,9000];
for m1 in [40,80,120,160,200,240,280,320,1000,2000,3000]
    global m = m1;
    global n,r = floor(Int64,0.5*m),floor(Int64,0.25*m);
    nameInst = string("A_",m,"_",n,"_",r);
    println(string("\nStarting instance: m,n,r: ",m,",",n,",",r))
    A = getMatlabInstance(nameInst,"A");
    inst = GinvInst(A,m,n,r);
    #ginvInit = getInitialInfoGinv(inst);
    U,S,V = svd(A,full=true)
    V1,V2 = V[:,1:r],V[:,r+1:n];
    U1,D = U[:,1:r],diagm(S[1:r]);
    V1Dinv = V1*inv(D);
    κ =  maximum(norm.(eachrow(V1)))
    Λ = 1/κ*V1
    # Λ = zeros(n,r)


    time_admm = @elapsed admmsol = runADMM21n(V1Dinv,V2,U1,Λ);
    admmsol.z = getnorm21(admmsol.H);
    admmsol.time = time_admm;
    admmres = getResultsOptimizer(inst,admmsol);
    writeCSV!(arr_resADMM,admmres,:ADMM_21)
    println(string("m = ",m,". Prob ADMM finished in ", admmsol.time, " sec with 2,1 norm ", admmsol.z));
    flush(stdout)

    # msksol = solve21msk(inst,ginvInit);
    # mskres = getResultsOptimizer(inst,msksol);
    # writeCSV!(arr_resMSK,mskres,:MSK)
    # println(string("m = ",m,". Prob MSK finished in ", msksol.time, " sec."));
    # flush(stdout)
    
    # grbsol = solve21grb(inst,ginvInit);
    # grbres = getResultsOptimizer(inst,grbsol);
    # writeCSV!(arr_resGRB,grbres,:GRB)
    # println(string("m = ",m,". Prob GRB finished in ", grbsol.time, " sec."));
    # flush(stdout)

    

    #writeAllCSV!(arr_res_all,mskres,grbres,admmres);
end

# n21E = 0.0
# E = rand(n,r)
# for i = (1:n)
#     d = zeros(n*r)
#     d[(r*(i-1)+1:r*i)] .= 1
#     D = diagm(d)
#     n21E += norm(D*vec(E'));
# end
# getnorm21(E),n21E



#getnorm21(H),min21norm(A)
# J = -V1*inv(D) + E - 1/mu *Lam

# model = Model(Gurobi.Optimizer)
# @variable(model,Z[1:n-r,1:r]);
# X = J - V2*Z
# @objective(model,Min, sum(X[i,j]^2 for i = (1:n), j = (1:r)))
# optimize!(model)
# Z = value.(Z)

