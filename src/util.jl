function getInitialInfoGinv(inst::GinvInst)
    U,Σ,V = svd(inst.A,full=true);
    D = Σ[1:inst.r];
    Dinv = 1 ./ D;
    U1,U2 = U[:,1:inst.r],U[:,inst.r+1:inst.m];
    V1,V2 = V[:,1:inst.r],V[:,inst.r+1:inst.n];
    V1Dinv = V1*SA.spdiagm(Dinv);
    V1DinvU1T = V1Dinv*U1';
    V2V2T = V2*V2';
    ginvInit = GinvInit(U1,U2,V1,V2,Dinv,V1Dinv,V1DinvU1T,V2V2T);
    return ginvInit
    
    
    
    
    # Θ21 = (1/maximum(norm.(eachrow(V1))))*V1
    # U1U1T,V2V2T = U1*U1',V2*V2'
    # input1 = GinvADMM1(U1,V2,Θ1,V1DinvU1T,U1U1T,V2V2T,TP)
    # input21 = GinvADMM21(U1,V2,Θ21,V1Dinv,V2V2T,TP)
    # input20 = GinvADMM20(U1,V2,V1Dinv,V2V2T,0,0.0)
    #return ginvInit
end

function getnorm0(A::Matrix{T},ϵ) where {T}
    norm_0 = 0.0;
    m,n = size(A);
    for i = (1:m)
        for j = (1:n)
            if abs(A[i,j]) > ϵ
                norm_0 += 1;
            end
        end
    end
    return norm_0
end

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

function getnorm21(A::Matrix{T}) where {T}
    m,n = size(A);
    norm_21 = sum(norm(A[i,:],2) for i = (1:m));
    return norm_21
end

function getnorm20(A::Matrix{T},ϵ) where {T}
    norm2rowsA = [norm(A[i, :]) for i in axes(A,1)]
    norm_20 = 0
    for el in norm2rowsA
        if el > ϵ
            norm_20 += 1
        end
    end
    return norm_20
end

function getZeroRows(A::Matrix{T},ϵ) where {T}
    m,n = size(A);
    count_rows0 = 0;
    for i = (1:m)
        flag = true;
        for j = (1:n)
            if (abs(A[i,j])) > ϵ
                flag = false;
                break
            end
        end
        if flag
            count_rows0 = count_rows0 + 1;
        end
    end
    return count_rows0
end



function getZeroCols(A::Matrix{T},ϵ) where {T}
    m,n = size(A);
    count_cols0 = 0;
    for j = (1:n)
        flag = true;
        for i = (1:m)
            if (abs(A[i,j])) > ϵ
                flag = false;
                break
            end
        end
        if flag
            count_cols0 = count_cols0 + 1;
        end
    end
    return count_cols0
end

    

function getResultsSolver(inst::GinvInst,sol::SolutionOptimizer)
    ϵ = 1e-5;
    A = inst.A;
    H = sol.H;
    norm_fro = norm(H,2);
    norm_21  = getnorm21(H);
    norm_0 = getnorm0(H,ϵ);
    norm_1 = getnorm1(H);
    # count_rows0 = getZeroRows(H,ϵ);
    # NZR = inst.n - count_rows0;
    NZR = getnorm20(H,ϵ)
    count_cols0 = getZeroCols(H,ϵ);
    NZC = inst.m - count_cols0;
    p1,p2,p3,p4 = norm(A*H*A - A),norm(H*A*H - H),norm(A*H - H'*A'),norm(H*A - A'*H');
    p1 < ϵ ? bool_p1 = true : bool_p1 = false;
    p2 < ϵ ? bool_p2 = true : bool_p2 = false;
    p3 < ϵ ? bool_p3 = true : bool_p3 = false;
    p4 < ϵ ? bool_p4 = true : bool_p4 = false;
    ginvResult = GinvResultSolver(
        inst.m,
        inst.n,
        inst.r,
        sol.z,
        NZC,
        NZR,
        norm_0,
        norm_1,
        norm_21,
        0, #iter
        sol.time,
        norm_fro,
        bool_p1,bool_p2,bool_p3,bool_p4,
        p1,p2,p3,p4
    )
    return ginvResult
end

function getResultsADMM(inst::GinvInst,sol::SolutionADMM)
    ϵ = 1e-5;
    A = inst.A;
    H = sol.H;
    norm_fro = norm(H,2);
    norm_21  = getnorm21(H);
    norm_0 = getnorm0(H,ϵ);
    norm_1 = getnorm1(H);
    # count_rows0 = getZeroRows(H,ϵ);
    # NZR = inst.n - count_rows0;
    NZR = getnorm20(H,ϵ);
    count_cols0 = getZeroCols(H,ϵ);
    NZC = inst.m - count_cols0;
    p1,p2,p3,p4 = norm(A*H*A - A),norm(H*A*H - H),norm(A*H - H'*A'),norm(H*A - A'*H');
    p1 < ϵ ? bool_p1 = true : bool_p1 = false;
    p2 < ϵ ? bool_p2 = true : bool_p2 = false;
    p3 < ϵ ? bool_p3 = true : bool_p3 = false;
    p4 < ϵ ? bool_p4 = true : bool_p4 = false;
    ginvResult = GinvResultADMM(
        inst.m,
        inst.n,
        inst.r,
        sol.z,
        NZC,
        NZR,
        norm_0,
        norm_1,
        norm_21,
        sol.iter,
        sol.time,
        norm_fro,
        sol.res_pri,sol.res_dual,sol.res_d_ML,sol.res_opt,sol.eps_p,sol.eps_d,
        bool_p1,bool_p2,bool_p3,bool_p4,
        p1,p2,p3,p4
    )
    return ginvResult
end

function getResultsExtendedADMM(inst::GinvInst,sol::SolutionExtendedADMM)
    ϵ = 1e-5;
    A = inst.A;
    H = sol.H;
    norm_fro = norm(H,2);
    norm_21  = getnorm21(H);
    norm_0 = getnorm0(H,ϵ);
    norm_1 = getnorm1(H);
    # count_rows0 = getZeroRows(H,ϵ);
    # NZR = inst.n - count_rows0;
    NZR = getnorm20(H,ϵ);
    count_cols0 = getZeroCols(H,ϵ);
    NZC = inst.m - count_cols0;
    p1,p2,p3,p4 = norm(A*H*A - A),norm(H*A*H - H),norm(A*H - H'*A'),norm(H*A - A'*H');
    p1 < ϵ ? bool_p1 = true : bool_p1 = false;
    p2 < ϵ ? bool_p2 = true : bool_p2 = false;
    p3 < ϵ ? bool_p3 = true : bool_p3 = false;
    p4 < ϵ ? bool_p4 = true : bool_p4 = false;
    ginvResult = GinvResultExtendedADMM(
        inst.m,
        inst.n,
        inst.r,
        sol.z,
        NZC,
        NZR,
        norm_0,
        norm_1,
        norm_21,
        sol.iter,
        sol.time,
        norm_fro,
        sol.res_cpE,sol.res_cdE,sol.res_cpB,sol.res_cdB,
        sol.eps_pE,sol.eps_dE,sol.eps_pB,sol.eps_dB,
        bool_p1,bool_p2,bool_p3,bool_p4,
        p1,p2,p3,p4
    )
    return ginvResult
end

function mkdir_ginv(folder_name)
    if !isdir(folder_name)
        mkdir(folder_name)
    end
end

function init_ginv_res_solver()
    ginv_res = [];
    ginv_label = [];
    for field in fieldnames(typeof(GinvResultSolver()))
        push!(ginv_label,string(field))
    end
    push!(ginv_res,ginv_label);
    return ginv_res
end
function init_ginv_res_admm()
    ginv_res = [];
    ginv_label = [];
    for field in fieldnames(typeof(GinvResultADMM()))
        push!(ginv_label,string(field))
    end
    push!(ginv_res,ginv_label);
    return ginv_res
end

function writeCSV!(arr_res::Vector{T},ginv_res::Union{GinvResultADMM,GinvResultExtendedADMM, GinvResultSolver},ginv_code::Symbol) where {T}
    folder_1 = "..\\ResultsCSV";
    mkdir_ginv(folder_1)
    # File name
    results_lieu   = string(folder_1,"\\",string("results_",ginv_code,".csv"));
    # Get infos
    new_ginv_vec = [];
    for field in fieldnames(typeof(ginv_res))
        push!(new_ginv_vec,getfield(ginv_res,field))
    end
    push!(arr_res,new_ginv_vec)
    # write information
    writedlm(results_lieu,arr_res, ',')
end

function round_exact(value,my_digits)
    val = string(round(value,digits = my_digits))
    count = 0;
    flag = false
    for el in val
        if !flag
            if el == '.'
                flag = true;
            end
        else
            count +=1;
        end
    end
    dif_float = (my_digits-count);
    return string(val,"0"^dif_float)
end

function append_to_plot2!(P1,D1,It,primal_res,dual_res,iter)
    primal_res <= 1e-15 ? primal_res = 1e-12 : primal_res
    dual_res <= 1e-15 ? dual_res = 1e-12 : dual_res
    push!(P1,primal_res)
    push!(D1,dual_res)
    push!(It,iter)
end

function append_to_plot_extended2!(primal_res1,primal_res2,dual_res1,dual_res2,iter)
    primal_res1 <= 1e-15 ? primal_res1 = 1e-12 : primal_res1
    primal_res2 <= 1e-15 ? primal_res2 = 1e-12 : primal_res2
    dual_res1 <= 1e-15 ? dual_res1 = 1e-12 : dual_res1
    dual_res2 <= 1e-15 ? dual_res2 = 1e-12 : dual_res2
    push!(P1,primal_res1)
    push!(P2,primal_res2)
    push!(D1,dual_res1)
    push!(D2,dual_res2)
    push!(It,iter)
end

# function get_input_21(ginvInit)
#     nput21 = GinvADMM21(U1,V2,Θ21,V1Dinv,V2V2T,TP)


# function writeAllCSV!(arr_res_all::Vector{T},ra::GinvResult,rb::GinvResult) where {T}
#     folder_1 = "..\\ResultsCSV";
#     mkdir_ginv(folder_1)
#     # File name
#     results_lieu   = string(folder_1,"\\",string("latex_results.csv"));
#     # Get infos
#     push!(arr_res_all,[ra.m,ra.n,ra.r,"A",ra.NZR,ra.norm_0,ra.norm_1,ra.norm_21,ra.iter,ra.time]);
#     push!(arr_res_all,[rb.m,rb.n,rb.r,"B",rb.NZR,rb.norm_0,rb.norm_1,rb.norm_21,rb.iter,rb.time]);
#     #push!(arr_res_all,[rc.m,rc.n,rc.r,"C",rc.NZR,rc.norm_0,rc.norm_1,rc.norm_21,rc.norm_fro,rc.time]);
#     # write information
#     writedlm(results_lieu,arr_res_all, ',')
# end

# function saveSolsJLD(inst::GinvInst,psols::ProbSols)
#     file_location = "..\\JLDFiles\\";
#     name_inst = string(
#         file_location,
#         inst.m,"_",
#         inst.n,"_",
#         inst.r,".jld"
#     );
#     JLD.save(name_inst, "psols",psols)
# end

