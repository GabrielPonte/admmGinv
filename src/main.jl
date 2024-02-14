using LinearAlgebra
using JuMP, MosekTools, Gurobi
using MAT,DelimitedFiles

import MutableArithmetics
const MA = MutableArithmetics
const SA = MA.SparseArrays

using JLD
using Plots, LaTeXStrings


include("types.jl")
include("util.jl")
include("instances.jl")
include("solver.jl")
include("admm.jl")

arr_resADMM_1,arr_resADMM_21,arr_resADMM_20,arr_resADMM_210 = init_ginv_res_admm(), init_ginv_res_admm(),init_ginv_res_admm(),init_ginv_res_admm();
arr_resGRB_1,arr_resMSK_21 = init_ginv_res_solver(),init_ginv_res_solver();
println("\nStarting procedure...")

TP = :ML
ω21 = 0.5
#arr_res_all = [];

global m,n,r  = 100,200,70;
M = [100,200,300,400,500,1000,2000,3000,4000,5000];

# for m1 in M
for m1 in [50]
    global m = m1;
    global n,r = floor(Int64,0.5*m),floor(Int64,0.25*m);
    nameInst = string("A_",m,"_",n,"_",r);
    println(string("\nStarting instance: m,n,r: ",m,",",n,",",r))
    # A = getMatlabInstance(nameInst,"A");
    global A = rand(m,r) * rand(r,n)
    global inst = GinvInst(A,m,n,r);
    # Initialization
    global U,S,V = svd(A,full=true)
    global V1,V2 = V[:,1:r],V[:,r+1:n];
    global U1,D = U[:,1:r],diagm(S[1:r]);
    global P1,D1,P2,D2,It = [],[],[],[],[];
    global V1Dinv,V1U1 = V1*inv(D),V1*U1';
    global G = V1Dinv*U1';
    global Λ1 = (1/maximum(abs.(V1U1)))*V1U1 
    global Λ21 = (1/maximum(norm.(eachrow(V1))))*V1
    
    time_admm_1 = @elapsed admmsol_1 = runADMM1n(G,V2,U1,Λ1,TP,false);
    admmsol_1.z = getnorm1(admmsol_1.H);
    admmsol_1.time = time_admm_1;
    admmres_1 = getResultsADMM(inst,admmsol_1);
    writeCSV!(arr_resADMM_1,admmres_1,Symbol(:ADMM_1_test,TP))
    println(string("m = ",m,". ADMM 1 finished in ", round_exact(admmsol_1.time,2), " sec. 1 norm ", round_exact(admmsol_1.z,3), ", 2,0 norm ", admmres_1.NZR, ". Iter: ", admmsol_1.iter));
    flush(stdout)
    
    time_admm_21 = @elapsed admmsol_21 = runADMM21n(V1Dinv,V2,U1,Λ21,TP,false);
    admmsol_21.z = getnorm21(admmsol_21.H);
    admmsol_21.time = time_admm_21;
    admmres_21 = getResultsADMM(inst,admmsol_21);
    writeCSV!(arr_resADMM_21,admmres_21,Symbol(:ADMM_21_,TP))
    println(string("m = ",m,". ADMM 2,1 finished in ", round_exact(admmsol_21.time,2), " sec. 2,1 norm ", round_exact(admmsol_21.z,3), ", 2,0 norm ", admmres_21.NZR, ". Iter: ", admmsol_21.iter));
    flush(stdout)


    
    # time_admm_20 = @elapsed admmsol_20 = runADMM20n(V1Dinv,V2,U1,admmres_21.NZR,ω21,false);
    # admmsol_20.z = getnorm21(admmsol_20.H);
    # admmsol_20.time = time_admm_20;
    # admmres_20 = getResultsADMM(inst,admmsol_20);
    # writeCSV!(arr_resADMM_20,admmres_20,Symbol(:ADMM_20,ω21))
    # println(string("m = ",m,". ADMM 2,0 finished in ", round_exact(admmsol_20.time,2), " sec. 2,1 norm ", round_exact(admmsol_20.z,3), ", 2,0 norm ", admmres_20.NZR,". Iter: ", admmsol_20.iter));
    # flush(stdout)

    # time_admm_210 = @elapsed admmsol_210 = runADMM210n(V1Dinv,V2,U1,Λ21,floor(0.75*n));
    # admmsol_210.z = getnorm21(admmsol_210.H);
    # admmsol_210.time = time_admm_210;
    # admmres_210 = getResultsExtendedADMM(inst,admmsol_210);
    # writeCSV!(arr_resADMM_210,admmres_210,:ADMM_210)
    # println(string("m = ",m,". ADMM 2,1,0 finished in ", round_exact(admmsol_210.time,2), " sec. 2,1 norm ", round_exact(admmsol_210.z,3), ", 2,0 norm ", admmres_210.NZR,". Iter: ", admmsol_210.iter));
    # flush(stdout)

    #plot(It,[P1,D1],labels=[L"\mathcal{P}" L"\mathcal{D}"], yaxis=:log, yticks =  10.0.^(5:-1:-15),xlabel = L"k" ,ylabel="residual")

    
    # if m <= 2000
    #     ginvInit = getInitialInfoGinv(inst);
    #     msksol = solve21msk(inst,ginvInit);
    #     mskres = getResultsSolver(inst,msksol);
    #     writeCSV!(arr_resMSK_21,mskres,:MSK_21)
    #     println(string("m = ",m,". Prob MSK 2,1 finished in ", msksol.time, " sec."));
    #     flush(stdout)
    # end
    # ginvInit = getInitialInfoGinv(inst);
    # grbsol = solve1grb(inst,ginvInit);
    # grbres = getResultsSolver(inst,grbsol);
    # writeCSV!(arr_resGRB_1,grbres,:GRB_1_TimeLimit)
    # println(string("m = ",m,". Prob GRB finished in ", grbsol.time, " sec."));
    # flush(stdout)

    

    #writeAllCSV!(arr_res_all,mskres,grbres,admmres);
end
# num_init = 1
# plot(It,[P1,D1],labels=[L"\mathcal{P}" L"\mathcal{D}"], yaxis=:log, yticks =  10.0.^(5:-1:-15),xlabel = L"k" ,ylabel="residual")
# plot(It,D1,label=string(ω21), xlabel =L"k",ylabel=L"{\|\|H\|\|}_{2,0}")

# plot(It,[P1,P2,D1,D2],labels=[L"\mathcal{P}_E" L"\mathcal{P}_B" L"\mathcal{D}_E" L"\mathcal{D}_B"], xlabel = L"k" ,ylabel="residual",yaxis=:log, yticks =  10.0.^(5:-1:-15))
