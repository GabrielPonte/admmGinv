using Pkg
Pkg.activate(".")

using LinearAlgebra
using JuMP, MosekTools, Gurobi
using DelimitedFiles, CSV, DataFrames
using SparseArrays

using JLD
using Plots, LaTeXStrings
using Printf

include("types.jl")
include("util.jl")
include("closed_form_sols.jl")
include("admm.jl")
include("solvers.jl")
include("local_search.jl")

arr_resADMM_1,arr_resADMM_0 = init_ginv_res_admm(),init_ginv_res_admm()
arr_resADMM_21,arr_resADMM_20,arr_resADMM_2120 = init_ginv_res_admm(),init_ginv_res_admm(),init_ginv_res_admm();
arr_resLSdet,arr_resLS21 = init_ginv_res_admm(),init_ginv_res_admm();
arr_resGRB_1,arr_resMSK_21 = init_ginv_res_solver(),init_ginv_res_solver();
println("\nStarting procedure...")


# Instances to test
instances_array = ["S1","S2","S3","S4","S5","L1","L2","L3","L4","L5",]; 

run_ls = false; # Run local search based on determinant
run_admm1n = false; # Run ADMM 1-norm for P1, P2, P3
run_admm21n = true;# Run ADMM 2,1-norm for P1, P2, P3
run_20n = true; # Run ADMM with target 2,0-norm
run_2120n = true; # Run ADMM-2,1 with target 2,0-norm
run_grb = false; # Run Gurobi solver for min 1-norm
run_msk = false; # Run Mosek solver for min 2,1-norm

fixed_tol = true # Choose between fixed tolerance for ADMM or dynamic tolerance proposed by Boyd

fixed_tol ? TP = :rel_tol : TP=:fixed_tol

for inst in instances_array
    A = Float64.(Matrix(CSV.read(string(".//Instances//",inst,".csv"), DataFrame; header=false)))
    m,n = size(A);
    r = floor(Int64,0.25*m);
    inst = GinvInst(A,m,n,r);
    if run_ls
        # Get linearly independent rows and columns
        R,C,time_RC = get_LI_RC(A,r)
        # Run local-search based on the determinant
        lsdetsol,C = LS_det(A,R,C,:FP)
        lsdetsol.time += time_RC 
        admmres_lsdet = getResultsADMM(inst,lsdetsol);
        writeCSV!(arr_resLSdet,admmres_lsdet,Symbol(:LS_det))
        @info (string("m = ",m,". LS det: ", round_exact(admmres_lsdet.time,2), " sec. 2,1 norm ", round_exact(admmres_lsdet.norm_21,3), ", 0 norm ", admmres_lsdet.norm_0, ". Iter: ", admmres_lsdet.iter));
        flush(stdout); GC.gc();
        # Run local-search based on the 2,1-norm
        ls21sol= LS_21(A,C,:FI);
        ls21sol.time += lsdetsol.time 
        admmres_ls21 = getResultsADMM(inst,ls21sol);
        writeCSV!(arr_resLS21,admmres_ls21,Symbol(:LS_21))
        @info (string("m = ",m,". LS 21: ", round_exact(admmres_ls21.time,2), " sec. 2,1 norm ", round_exact(admmres_ls21.norm_21,3), ", 0 norm ", admmres_ls21.norm_0, ". Iter: ", admmres_ls21.iter));
        flush(stdout); GC.gc();
    end
    # ADMM Initialization
    ginvInit = getInitialInfoGinv(inst)
    if run_admm1n
        if fixed_tol
            # Run ADMM based on the 1-norm
            time_admm_1 = @elapsed admmsol_1 = admm1norm(ginvInit;eps_opt=1e-2,stop_limit=:OptGap);
        else
            # Run ADMM based on the 1-norm
            time_admm_1 = @elapsed admmsol_1 = admm1norm(ginvInit;eps_abs=1e-4,eps_rel=1e-4,stop_limit=:Boyd);
        end
        admmres_1 = getResultsADMM(inst,admmsol_1);
        writeCSV!(arr_resADMM_1,admmres_1,Symbol(:ADMM_1_,TP))
        @info (string("m = ",m,". ADMM 1: ", round_exact(admmsol_1.time,2), " sec. 1 norm ", round_exact(admmsol_1.z,3), ", 0 norm ", admmres_1.norm_0, ". Iter: ", admmsol_1.iter));
        flush(stdout); GC.gc();
    end
    # Run ADMM based on the 2,1-norm
    if run_admm21n || run_2120n || run_20n
        time_admm_21 = @elapsed admmsol_21,M_21,Λ_21 = admm21norm(ginvInit);
        admmres_21 = getResultsADMM(inst,admmsol_21);
        if run_admm21n
            writeCSV!(arr_resADMM_21,admmres_21,Symbol(:ADMM_21_,TP))
            @info (string("m = ",m,". ADMM 2,1: ", round_exact(admmsol_21.time,2), " sec. |2,1|: ", round_exact(admmsol_21.z,3), ". |2,0|: ", admmres_21.norm_20, ". |0|: ", admmres_21.norm_0,". Iter: ", admmsol_21.iter));
            flush(stdout); GC.gc();
        end
        if run_2120n || run_20n
            # Run ADMM based on the 2,1-2,0-norm
            if m >= 1000
                for ω21 in [0.25, 0.50, 0.75, 0.80, 0.90, 0.95]
                    nzr21 = admmres_21.norm_20;
                    if run_20n
                        time_admm_20 = @elapsed admmsol_20 = admm20norm(ginvInit,ω21,nzr21);
                        admmres_20 = getResultsADMM(inst,admmsol_20);
                        writeCSV!(arr_resADMM_20,admmres_20,Symbol(:ADMM_20_,trunc(Int64,100*ω21)))
                        @info (string("w = ",round_exact(ω21,2),". ADMM 2,0: ", round_exact(admmsol_20.time,2), " sec. |2,1|: ", round_exact(admmsol_20.z,3), ". |2,0|: ", admmres_20.norm_20,". |0|: ", admmres_20.norm_0,". Iter: ", admmsol_20.iter));
                        flush(stdout); GC.gc();
                    end
                    if run_2120n
                        time_admm_2120 = @elapsed admmsol_2120 = admm2120norm(ginvInit,ω21,nzr21,M_21,Λ_21);
                        admmres_2120 = getResultsADMM(inst,admmsol_20);
                        writeCSV!(arr_resADMM_2120,admmres_2120,Symbol(:ADMM_2120_,trunc(Int64,100*ω21)))
                        @info (string("w = ",round_exact(ω21,2),". ADMM 2,0: ", round_exact(admmsol_2120.time,2), " sec. |2,1|: ", round_exact(admmsol_2120.z,3), ". |2,0|: ", admmres_2120.norm_20,". |0|: ", admmres_2120.norm_0,". Iter: ", admmsol_2120.iter));
                        flush(stdout); GC.gc();
                    end
                end
            end
        end

    end
    if run_grb 
        time_grb_1 = @elapsed grbsol_1 =solve1grb(inst,ginvInit)
        grbres_1 = getResultsSolver(inst,grbsol_1)
        writeCSV!(arr_resGRB_1,grbres_1,Symbol(:GRB_1))
        @info (string("m = ",m,". GRB 1: ",  round_exact(grbres_1.time,2), " sec. 1 norm ", round_exact(grbres_1.norm_1,3), ", 0 norm ", grbres_1.norm_0));
        flush(stdout); GC.gc();
    end
    if run_msk
        time_msk_21 = @elapsed msksol_21 =solve21msk(inst,ginvInit)
        mskres_21 = getResultsSolver(inst,msksol_21)
        writeCSV!(arr_resMSK_21,mskres_21,Symbol(:MSK_21))
        @info (string("m = ",m,". MSK 2,1: ", round_exact(mskres_21.time,2), " sec. |2,1|: ", round_exact(mskres_21.norm_21,3), ". |2,0|: ", mskres_21.norm_20, ". |0|: ", mskres_21.norm_0));
        flush(stdout); GC.gc();
    end
end