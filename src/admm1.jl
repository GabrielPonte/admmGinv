using LinearAlgebra, JuMP,MosekTools, Gurobi
using LinearAlgebra
using MAT,DelimitedFiles

import MutableArithmetics
const MA = MutableArithmetics
const SA = MA.SparseArrays

using JLD
#using Profile

include("types.jl")
include("util.jl")
include("probs.jl")
include("instances.jl")


function getnorm1(A::Matrix{T}) where {T}
    norm_1 = 0.0;
    m,n = size(A);
    for i = (1:m)
        for j = (1:n)
            norm_1 += abs(A[i,j]);
        end
    end
    return norm_1
end

function solve1grb(inst::GinvInst,ginvInit::GinvInit)
    m,n,r = inst.m,inst.n,inst.r;
    model = Model(Gurobi.Optimizer);
    set_silent(model);
    set_time_limit_sec(model, 7200)
    # set_attribute(model, "BarConvTol", 1e-5)
    # set_attribute(model, "BarQCPConvTol", 1e-5)
    # set_attribute(model, "FeasibilityTol", 1e-5)
    # set_attribute(model, "OptimalityTol", 1e-5)
    @variable(model, F[1:n,1:m]);
    @variable(model, Z[1:n-r,1:r]);
    Amat = ginvInit.G + ginvInit.Vs*Z*ginvInit.Ur';
    @constraint(model, F .>= Amat);
    @constraint(model, F .>= -Amat);
    @objective(model,Min,sum(F))
    time_to_solve = @elapsed optimize!(model);
    bsol = SolutionOptimizer();
    primal_status = termination_status(model);
    if primal_status != MOI.TIME_LIMIT
        bsol.H = ginvInit.G + ginvInit.Vs*value.(Z)*ginvInit.Ur';
        bsol.z,bsol.time =  objective_value(model),time_to_solve;
        bsol.iter = 0;
    else
        bsol.H = zeros(n,m);
        bsol.z,bsol.time = Inf, time_to_solve;
        bsol.iter = 0;
    end
    return bsol
end



function update_E1n(G,W,mu,Λ)
    mu_inv = 1/mu
    Y = G + W + mu_inv*Λ
    # Y = G + W + mu_inv*Λ
    E = sign.(Y).*max.(abs.(Y) .- (mu_inv), 0)
    # E = zeros(n,m)
    # for i in (1:n)
    #     for j in (1:m)
    #         y_ij = Y[i,j]
    #         if abs(y_ij) > mu_inv
    #             if y_ij > mu_inv
    #                 E[i,j] = y_ij - mu_inv;
    #             elseif y_ij < -mu_inv
    #                 E[i,j] = y_ij + mu_inv
    #             end
    #         end
    #     end
    # end
    return E
end

function update_Z1n(G,V2,U1,E,mu,Λ)
    return V2'*(- G + E - (1/mu)*Λ)*U1
end

function runADMM1n(G,V2,U1,Λ)
    mu = 1e-3
    # Λ = zeros(n,m)
    # E = zeros(n,r)
    # Z = zeros(n-r,r)
    W = zeros(n,m)
    # 
    V2V2T = V2*V2';
    U1U1T = U1*U1';
    iter = 0
    ρ = 1.01
    rho_iter = ρ
    pri_updt,dual_updt,usual_updt = 0,0,0
    while true
        E = update_E1n(G,W,mu,Λ)
        J = (- G + E - (1/mu)*Λ);
        # Z = update_Z1n(G,V2,U1,E,mu,Λ)
        W = V2V2T*J*U1U1T 
        res_infeas = G + W - E
        primal_res = norm(res_infeas)
        # dual_res = maximum([norm(V2'*Λ*U1),maximum(abs.(Λ))-1 ])
        dual_res = maximum(abs.(Λ))-1
        Λ = (Λ + mu*res_infeas)
        # mu = ρ*mu
        #@show dual_res
        # @show norm(primal_res), dual_res
        if primal_res > 10*dual_res
            #@show "oi"
            mu = 1.1*mu 
            # rho_iter = 1.1
            pri_updt += 1
        elseif dual_res > 10*primal_res
            #@show maximum(abs.(Λ)),maximum(abs.(Λ))-1
            mu = 1.005*mu
            # rho_iter = 1.005
            dual_updt += 1
        else
            mu = ρ*mu
            # rho_iter = ρ
            usual_updt += 1
        end
        iter += 1
        
        #@show maximum(abs.(Λ))
        # @show primal_res,norm(V2'*Λ*U1),maximum(abs.(Λ))
        # @show mu
        if (primal_res < 1e-5 && dual_res <= 1.01) || iter == 5000
            if !(norm(V2'*Λ*U1) >= 1e-5 && iter < 5000)
                # @show primal_res
                # @show norm(V2'*Λ*U1) 
                # @show maximum(abs.(Λ))
                # @show pri_updt,dual_updt,usual_updt
                break
            end
        end
    end
    println(string("Number of iterations: ", iter))
    # @show Z
    admmsol = SolutionOptimizer();
    admmsol.H = G + W;
    admmsol.iter = iter;
    return admmsol
end

arr_res_all = [];
arr_resADMM,arr_resGRB,arr_resMSK = init_ginv_res(),init_ginv_res(),init_ginv_res();
dens = 1.00
global m,n,r  = 100,200,70;
# M = [40,80,120,160,200,240,280,320,1000,2000,3000,4000,5000,6000,7000,8000,9000];
M = [40,80,120,160,200,240,280,320,1000,2000,3000,4000,5000];
for m1 in M
    global m = m1;
    global n,r = floor(Int64,0.5*m),floor(Int64,0.25*m);
    nameInst = string("A_",m,"_",n,"_",r);
    println(string("\nStarting instance: m,n,r: ",m,",",n,",",r))
    A = getMatlabInstance(nameInst,"A");
    inst = GinvInst(A,m,n,r);
    ginvInit = getInitialInfoGinv(inst);
    U,S,V = svd(A,full=true)
    V1 = V[:,1:r];
    V2 = V[:,r+1:n];
    U1 = U[:,1:r];
    global S = S;
    D = diagm(S[1:r]);
    V1Dinv = V1*inv(D);
    G = V1Dinv*U1';
    # Z = -V2'*V1*inv(D) 
    V1U1 = V1*U1'
    Λ = (1/maximum(V1U1))*V1U1 
    # my_time = @elapsed @profview runADMM1n(G,V2,U1);
    # @show my_time
    time_admm = @elapsed admmsol = runADMM1n(G,V2,U1,Λ);
    admmsol.z = getnorm1(admmsol.H);
    admmsol.time = time_admm;
    admmres = getResultsOptimizer(inst,admmsol);
    writeCSV!(arr_resADMM,admmres,:ADMM_1)
    println(string("m = ",m,". Prob ADMM finished in ", admmsol.time, " sec with 1 norm ", admmsol.z));
    #println(admmsol.iter)
    flush(stdout)

    # msksol = solve21msk(inst,ginvInit);
    # mskres = getResultsOptimizer(inst,msksol);
    # writeCSV!(arr_resMSK,mskres,:MSK)
    # println(string("m = ",m,". Prob MSK finished in ", msksol.time, " sec."));
    # flush(stdout)
    
    # grbsol = solve1grb(inst,ginvInit);
    # grbres = getResultsOptimizer(inst,grbsol);
    # writeCSV!(arr_resGRB,grbres,:GRB_1)
    # println(string("m = ",m,". Prob GRB finished in ", grbsol.time, " sec."));
    # flush(stdout)

    

    # writeAllCSV!(arr_res_all,grbres,admmres);
end





#getnorm21(H),min21norm(A)
# J = -V1*inv(D) + E - 1/mu *Lam

# model = Model(Gurobi.Optimizer)
# @variable(model,Z[1:n-r,1:r]);
# X = J - V2*Z
# @objective(model,Min, sum(X[i,j]^2 for i = (1:n), j = (1:r)))
# optimize!(model)
# Z = value.(Z)

