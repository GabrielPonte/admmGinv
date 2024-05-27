using LinearAlgebra
using JuMP, MosekTools, Gurobi
using MAT,DelimitedFiles

import MutableArithmetics
const MA = MutableArithmetics
const SA = MA.SparseArrays
using Distributed

using JLD
using Plots, LaTeXStrings

using Profile

include("types.jl")
include("util.jl")
include("instances.jl")
include("closed_form_sols.jl")
include("admm.jl")
include("step_sizes.jl")
# include("solvers.jl")
include("local_search.jl")

arr_resADMM_1,arr_resADMM_0 = init_ginv_res_admm(),init_ginv_res_admm()
arr_resADMM_21,arr_resADMM_20,arr_resADMM_210 = init_ginv_res_admm(),init_ginv_res_admm(),init_ginv_res_admm();
arr_resBREG_1 = init_ginv_res_admm();
arr_resGRB_1,arr_resMSK_21 = init_ginv_res_solver(),init_ginv_res_solver();
println("\nStarting procedure...")

TP = :ML
ω21 = 0.5
#arr_res_all = [];




global m,n,r  = 100,200,70;
M = [100,200,300,400,500,1000,2000,3000,4000,5000];

pres,dres,tols,objs,rhos = [],[],[],[],[];

# for m1 in [200]#[1000,2000,3000,4000,5000]
m1 = 1000;
m = m1;
n,r = floor(Int64,0.5*m),floor(Int64,0.25*m);
nameInst = string("A_",m,"_",n,"_",r);
println(string("\nStarting instance: m,n,r: ",m,",",n,",",r))
A = getMatlabInstance(nameInst,"A");
R,C,time_RC = get_LI_RC(A,r)
H,time_det,swaps_det,det_Ar,C = LS_det(A,R,C,:FP)
n21old = round(getnorm21(H),digits=3)
@show time_det
# BLAS.set_num_threads(1) # 1 is better than 32 (?)
# println("32")
@profview LS_21(A,C,:FI);
# H,time_ls21,swaps = LS_21(A,C,:FI);
# swaps
# n21new = round(getnorm21(H),digits=3)
# n21old,n21new
    # @show swaps_det,swaps_21
    # @show rank(A[R,C],1e-6),r
    # global A = rand(m,r) * rand(r,n)
    # inst = GinvInst(A,m,n,r);
    # Initialization
    # ginvInit = getInitialInfoGinv(inst)

    # time_admm_1 = @elapsed admmsol_1 = admm1norm_v2(ginvInit,eps_opt=1e-4,stop_limit=:OptGap);
    # admmsol_1.z = getnorm1(admmsol_1.H);
    # admmsol_1.time = time_admm_1;
    # admmres_1 = getResultsADMM(inst,admmsol_1);
    # writeCSV!(arr_resADMM_1,admmres_1,Symbol(:ADMM_1_test,TP))
    # @show admmsol_1.z,admmsol_1.time
    # println(string("m = ",m,". ADMM 1 finished in ", round_exact(admmsol_1.time,2), " sec. 1 norm ", round_exact(admmsol_1.z,3), ", 0 norm ", admmres_1.norm_0, ". Iter: ", admmsol_1.iter));
    # flush(stdout)

    # time_breg_1 = @elapsed bregsol_1 = breg1norm(ginvInit,eps=5e-1);
    # bregsol_1.z = getnorm1(bregsol_1.H);
    # bregsol_1.time = time_breg_1;
    # bregres_1 = getResultsADMM(inst,bregsol_1);
    # writeCSV!(arr_resBREG_1,bregres_1,Symbol(:BREG_1,TP))
    # println(string("m = ",m,". BREG 1 finished in ", round_exact(bregsol_1.time,2), " sec. 1 norm ", round_exact(bregsol_1.z,3), ", 0 norm ", bregres_1.norm_0, ". Iter: ", bregsol_1.iter));
    # flush(stdout)
    
    # time_admm_21 = @elapsed admmsol_21,M_21,Λ_21 = admm21norm(ginvInit);#runADMM21n(V1Dinv,V2,U1,Λ21,TP,false);
    # admmsol_21.z = getnorm21(admmsol_21.H);
    # admmsol_21.time = time_admm_21;
    # admmres_21 = getResultsADMM(inst,admmsol_21);
    # writeCSV!(arr_resADMM_21,admmres_21,Symbol(:ADMM_21_,TP))
    # println(string("m = ",m,". ADMM 2,1: ", round_exact(admmsol_21.time,2), " sec. |2,1|: ", round_exact(admmsol_21.z,3), ". |2,0|: ", admmres_21.NZR, ". |0|: ", admmres_21.norm_0,". Iter: ", admmsol_21.iter));
    # flush(stdout)

    # time_admm_0 = @elapsed admmsol_0 = admm0norm(ginvInit,trunc(Int64,m*(n-2)));
    # admmsol_0.z = getnorm1(admmsol_0.H);
    # admmsol_0.time = time_admm_0;
    # admmres_0 = getResultsADMM(inst,admmsol_0);
    # writeCSV!(arr_resADMM_0,admmres_0,Symbol(:ADMM_0_test,TP))
    # @show admmsol_0.z,admmsol_0.time
    # println(string("m = ",m,". ADMM 0 finished in ", round_exact(admmsol_0.time,2), " sec. 1 norm ", round_exact(admmsol_0.z,3), ", 0 norm ", admmres_0.norm_0, ". Iter: ", admmsol_0.iter));
    # flush(stdout)

    # """
    # 0.70 = 1e2 
    # 0.90 = 1e3
    # 0.95 = 1e4

    # """
    # ω21 = 0.95; nzr21 = admmres_21.NZR;
    # time_admm_20 = @elapsed admmsol_20 = admm20norm(ginvInit,ω21,nzr21,M_21,Λ_21;rho=1e3);#runADMM20n(V1Dinv,V2,U1,admmres_21.NZR,ω21,false);
    # admmsol_20.z = getnorm21(admmsol_20.H);
    # admmsol_20.time = time_admm_20;
    # admmres_20 = getResultsADMM(inst,admmsol_20);
    # writeCSV!(arr_resADMM_20,admmres_20,Symbol(:ADMM_2120_,trunc(Int64,100*ω21)))
    # println(string("m = ",m,". ADMM 2,0: ", round_exact(admmsol_20.time,2), " sec. |2,1|: ", round_exact(admmsol_20.z,3), ". |2,0|: ", admmres_20.NZR,". |0|: ", admmres_20.norm_0,". Iter: ", admmsol_20.iter));
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
# end

# num_init = 1
# plot(It,[P1,D1],labels=[L"\mathcal{P}" L"\mathcal{D}"], yaxis=:log, yticks =  10.0.^(5:-1:-15),xlabel = L"k" ,ylabel="residual")
# plot(It,D1,label=string(ω21), xlabel =L"k",ylabel=L"{\|\|H\|\|}_{2,0}")

# plot(It,[P1,P2,D1,D2],labels=[L"\mathcal{P}_E" L"\mathcal{P}_B" L"\mathcal{D}_E" L"\mathcal{D}_B"], xlabel = L"k" ,ylabel="residual",yaxis=:log, yticks =  10.0.^(5:-1:-15))

# y_vals = rhos
# # y_vals = pres
# # y_vals = dres
# x_vals = collect(1:length(y_vals))

# plot(
#         x_vals,
#         y_vals,
#         lw = 1.25,
#         yaxis=:log,
#         legend=:topright,
#         # legendfontsize=4,
#         # palette=palette(:default)#[init_color:end],
#         #palette = palette(:darkrainbow)[[1,3]]
#     )
