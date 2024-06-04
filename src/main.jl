using LinearAlgebra
using JuMP, MosekTools, Gurobi
using MAT,DelimitedFiles

import MutableArithmetics
const MA = MutableArithmetics
const SA = MA.SparseArrays
# using Base.Threads

using JLD
using Plots, LaTeXStrings
using Printf

include("types.jl")
include("util.jl")
include("closed_form_sols.jl")
include("admm.jl")
# include("solvers.jl")
include("local_search.jl")

arr_resADMM_1,arr_resADMM_0 = init_ginv_res_admm(),init_ginv_res_admm()
arr_resADMM_21,arr_resADMM_20,arr_resADMM_210 = init_ginv_res_admm(),init_ginv_res_admm(),init_ginv_res_admm();
arr_resLSdet,arr_resLS21 = init_ginv_res_admm(),init_ginv_res_admm();
arr_resGRB_1,arr_resMSK_21 = init_ginv_res_solver(),init_ginv_res_solver();
println("\nStarting procedure...")

TP = :ML

global m,n,r  = 100,200,70;
M2 = [100,200,300,400,500,1000,2000,3000,4000,5000];
dens = 1.00;

for m1 in [1000]#M2
    m = m1;
    n,r = floor(Int64,0.5*m),floor(Int64,0.25*m);
    dens_int = @sprintf("%03d", trunc(Int64,dens*100));
    nameInst = string("A",dens_int,"_",m,"_",n,"_",r);
    println(string("\nStarting instance: m,n,r: ",m,",",n,",",r))
    A = getMatlabInstance(nameInst,"A");
    inst = GinvInst(A,m,n,r);
    # Get linearly independent rows and columns
    R,C,time_RC = get_LI_RC(A,r)
    # Run local-search based on the determinant
    lsdetsol,C = LS_det(A,R,C,:FP)
    lsdetsol.time += time_RC 
    admmres_lsdet = getResultsADMM(inst,lsdetsol);
    println(string("m = ",m,". LS det: ", round_exact(admmres_lsdet.time,2), " sec. 2,1 norm ", round_exact(admmres_lsdet.norm_21,3), ", 0 norm ", admmres_lsdet.norm_0, ". Iter: ", admmres_lsdet.iter));
    writeCSV!(arr_resLSdet,admmres_lsdet,Symbol(:LS_det_,dens_int))
    # Run local-search based on the 2,1-norm
    ls21sol= LS_21(A,C,:FI);
    ls21sol.time += lsdetsol.time 
    admmres_ls21 = getResultsADMM(inst,ls21sol);
    println(string("m = ",m,". LS 21: ", round_exact(admmres_ls21.time,2), " sec. 2,1 norm ", round_exact(admmres_ls21.norm_21,3), ", 0 norm ", admmres_ls21.norm_0, ". Iter: ", admmres_ls21.iter));
    writeCSV!(arr_resLS21,admmres_ls21,Symbol(:LS_21_,dens_int))
    # ADMM Initialization
    ginvInit = getInitialInfoGinv(inst)
    # Run ADMM based on the 1-norm
    time_admm_1 = @elapsed admmsol_1 = admm1norm(ginvInit);
    admmres_1 = getResultsADMM(inst,admmsol_1);
    writeCSV!(arr_resADMM_1,admmres_1,Symbol(:ADMM_1_,TP,:_,dens_int))
    println(string("m = ",m,". ADMM 1: ", round_exact(admmsol_1.time,2), " sec. 1 norm ", round_exact(admmsol_1.z,3), ", 0 norm ", admmres_1.norm_0, ". Iter: ", admmsol_1.iter));
    flush(stdout)
    # Run ADMM based on the 2,1-norm
    time_admm_21 = @elapsed admmsol_21,M_21,Λ_21 = admm21norm(ginvInit);
    admmres_21 = getResultsADMM(inst,admmsol_21);
    writeCSV!(arr_resADMM_21,admmres_21,Symbol(:ADMM_21_,TP,:_,dens_int))
    println(string("m = ",m,". ADMM 2,1: ", round_exact(admmsol_21.time,2), " sec. |2,1|: ", round_exact(admmsol_21.z,3), ". |2,0|: ", admmres_21.norm_20, ". |0|: ", admmres_21.norm_0,". Iter: ", admmsol_21.iter));
    flush(stdout)
    # Run ADMM based on the 2,1-2,0-norm
    if m >= 1000
        for ω21 in [0.25, 0.50, 0.75, 0.80, 0.90, 0.95]
            nzr21 = admmres_21.norm_20;
            time_admm_20 = @elapsed admmsol_20 = admm2120norm(ginvInit,ω21,nzr21,M_21,Λ_21);
            admmres_20 = getResultsADMM(inst,admmsol_20);
            writeCSV!(arr_resADMM_20,admmres_20,Symbol(:ADMM_2120_,trunc(Int64,100*ω21),:_,dens_int))
            println(string("w = ",round_exact(ω21,2),". ADMM 2,0: ", round_exact(admmsol_20.time,2), " sec. |2,1|: ", round_exact(admmsol_20.z,3), ". |2,0|: ", admmres_20.norm_20,". |0|: ", admmres_20.norm_0,". Iter: ", admmsol_20.iter));
            flush(stdout)
        end
    end
end