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
include("closed_form_sols.jl")
include("admm.jl")
include("step_sizes.jl")
include("solvers.jl")

arr_resADMM_1,arr_resADMM_21,arr_resADMM_20,arr_resADMM_210 = init_ginv_res_admm(), init_ginv_res_admm(),init_ginv_res_admm(),init_ginv_res_admm();
arr_resGRB_1,arr_resMSK_21 = init_ginv_res_solver(),init_ginv_res_solver();
println("\nStarting procedure...")

TP = :ML
ω21 = 0.5
#arr_res_all = [];

global m,n,r  = 100,200,70;
# M = [100,200,300,400,500,1000,2000,3000,4000,5000];
M = [100,200,300,400,500,1000,2000,3000,4000,5000];


pres,dres,tols,objs,rhos = [],[],[],[],[];
for m1 in [40,100,200]
# for m1 in [1000]
    m = m1;
    n,r = floor(Int64,0.5*m),floor(Int64,0.25*m);
    nameInst = string("A_",m,"_",n,"_",r);
    println(string("\nStarting instance: m,n,r: ",m,",",n,",",r))
    A = getMatlabInstance(nameInst,"A");
    # global A = rand(m,r) * rand(r,n)
    inst = GinvInst(A,m,n,r);
    # Initialization
    ginvInit = getInitialInfoGinv(inst)
    time_grb_1 = @elapsed grbsol_1 = solve1_not_eff_grb(inst)#admm1norm(ginvInit);
    grbsol_1.z = getnorm1(grbsol_1.H);
    grbsol_1.time = time_grb_1;
    grbres_1 = getResultsSolver(inst,grbsol_1);
    writeCSV!(arr_resADMM_1,grbres_1,Symbol(:GRB_2_not_eff))
    @show grbsol_1.z,grbsol_1.time
    println(string("m = ",m,". ADMM 1 finished in ", round_exact(grbsol_1.time,2), " sec. 1 norm ", round_exact(grbsol_1.z,3)));
    flush(stdout)

    # bsol = solve1grb(inst,ginvInit)
    # @show bsol.z
    
    if time_grb_1 >= 7200
        break
    end


    

    #writeAllCSV!(arr_res_all,mskres,grbres,admmres);
end
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