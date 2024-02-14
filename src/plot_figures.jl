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

TP = :OTM

function get_plot(x_vals,y_vals,my_label,x_label,y_label)
    tot_procedures = 5;
    tam_label = length(my_label);
    init_color = (1 + tot_procedures - tam_label);
    return plot(
        x_vals,
        y_vals,
        labels = my_label;
        lw = 1.25,
        xlabel = x_label,
        # yaxis=:log,
        #xticks = 50:25:200,
        ylabel = y_label,
        # yticks =  10.0.^(5:-1:-15),
        #marker=(:circle,2),
        legend=:topright,
        # legendfontsize=4,
        # palette=palette(:default)#[init_color:end],
        palette = palette(:darkrainbow)[[1,3]]
    )
end

function get_plot20(x_vals,y_vals,my_label,x_label,y_label)
    tot_procedures = 5;
    tam_label = length(my_label);
    init_color = (1 + tot_procedures - tam_label);
    return plot(
        x_vals,
        y_vals,
        labels = my_label;
        lw = 1,
        xlabel = x_label,
        #yaxis=:log,
        #xticks = 50:25:200,
        ylabel = y_label,
        #marker=(:circle,2),
        legend=:topright,
        palette=palette(:default)#palette(:darktest)#[1:end],
    )
end

function save_plot20(x_vals,y_vals,my_label,x_label,y_label,name)
    p1 = get_plot20(x_vals,y_vals,my_label,x_label,y_label)
    folder_name = "Plots"
    if !isdir(string(".\\",folder_name))
        mkdir(string(".\\",folder_name))
    end
    savefig(p1,string(".\\Plots\\",name,".pdf"))
end

function save_plot(x_vals,y_vals,my_label,x_label,y_label,name)
    p1 = get_plot(x_vals,y_vals,my_label,x_label,y_label)
    folder_name = "Plots"
    if !isdir(string(".\\",folder_name))
        mkdir(string(".\\",folder_name))
    end
    savefig(p1,string(".\\Plots\\",name,".pdf"))
end

function runADMM1nplot(G,V2,U1,Θ,TP::Symbol)
    # initial params
    max_iter = 5e4;
    eps_abs,eps_rel = 1e-4,1e-4;
    rho = 2;
    rho_inv = 1/rho
    V2V2T = V2*V2';
    U1U1T = U1*U1';
    iter = 0;
    primal_res,dual_res = 0.0,0.0;
    Λ = rho_inv*Θ;
    E = G+Λ;
    E_old = E;    
    GW = Nothing
    if TP == :ML
        eps_p,eps_d = 0.0,0.0
    else
        eps_p,eps_d = 1e-4,1e-4*rho_inv
    end
    while true
        # update Z where W := V2*Z*U1'
        J = (- G + E - Λ)
        W = V2V2T*J*U1U1T
        # Update E
        GW = G + W
        Y = GW + Λ
        E = sign.(Y).*max.(abs.(Y) .- (rho_inv), 0)
        # Compute residuals
        res_infeas = GW - E
        primal_res = norm(res_infeas)
        # Update P
        Λ += res_infeas
        # Stopping criteria
        iter += 1
        if TP == :ML
            eps_p = sqrt(n*m)*eps_abs + eps_rel*maximum([norm(E),norm(GW)])
            eps_d = sqrt(n*m)*eps_abs + eps_rel*rho*norm(Λ)
            dual_res = rho*norm(E-E_old)
        else
            dual_res = norm(V2'*Λ*U1)
        end
        # plots
        append_to_plot!(P1,D1,It,primal_res,dual_res,iter)
        if (primal_res <= eps_p && dual_res <= eps_d) || iter == max_iter
            # @show primal_res,dual_res,iter
            break
        end
        E_old = E
    end
    Θ = rho*Λ
    admmsol = SolutionADMM();
    admmsol.H = GW;
    admmsol.iter = iter;
    admmsol.res_cp,admmsol.res_cd1,admmsol.res_cd2 = primal_res,dual_res,norm(V2'*Θ*U1);
    admmsol.eps_p,admmsol.eps_d = eps_p,eps_d/rho; 
    return admmsol
end

function runADMM21nplot(V1Dinv,V2,U1,Θ,TP::Symbol)
    # initial params
    max_iter = 5000;
    rho = 1;
    rho_inv = 1/rho
    iter = 0;
    primal_res,dual_res = 0.0,0.0;
    V2V2T = V2*V2'
    Λ = rho_inv*Θ;
    E = V1Dinv + Λ;
    E_old = E;  
    GW = Nothing  
    if TP == :ML
        eps_abs,eps_rel = 1e-6,1e-6;
        eps_p,eps_d = 0.0,0.0
    else
        eps_p,eps_d = 1e-6,1e-6*rho_inv
    end
    while true
        V1DinvΛ = V1Dinv + Λ
        # update Z where W := V2*Z
        W  = V2V2T*(-V1DinvΛ + E)
        # update E
        Y = V1DinvΛ + W
        E = zeros(n,r)
        norm_rows_y = norm.(eachrow(Y))
        for (i,el) in enumerate(norm_rows_y)
            if rho_inv < el
                E[i,:] = ((el - rho_inv)/el) * Y[i,:]
            end
        end
        # Compute residuals
        GW = V1Dinv + W
        res_infeas = GW -E
        primal_res = norm(res_infeas)
        # Update P
        Λ += res_infeas        
        # Stopping criteria
        iter += 1
        if TP == :ML
            eps_p = sqrt(n*r)*eps_abs + eps_rel*maximum([norm(E),norm(GW)])
            eps_d = sqrt(n*r)*eps_abs + eps_rel*rho*norm(Λ)
            dual_res = rho*norm(E-E_old)
        else
            dual_res = norm(V2'*Λ)
        end
        # plots
        append_to_plot!(P1,D1,It,primal_res,dual_res,iter)
        if (primal_res <= eps_p && dual_res <= eps_d) || iter == max_iter
            # @show primal_res,dual_res,iter
            break
        end
        E_old = E
    end
    Θ = rho*Λ
    admmsol = SolutionADMM();
    admmsol.H = GW*U1';
    admmsol.iter = iter;
    admmsol.res_cp,admmsol.res_cd1,admmsol.res_cd2 = primal_res,dual_res,norm(V2'*Θ);
    admmsol.eps_p,admmsol.eps_d = eps_p,eps_d/rho; 
    return admmsol
end

function runADMM20nplot(V1Dinv,V2,U1,nzr21,w21)
    # initial params
    max_iter = 5000;
    
    rho = 1;
    rho_inv = 1/rho
    eps_abs,eps_rel = 1e-7,1e-7*rho_inv;
    V2V2T = V2*V2';
    iter = 0;
    primal_resB,dual_resB = 0.0,0.0;
    eps_pB,eps_dB = 0.0,0.0;
    Φ,B = zeros(n,r),V1Dinv
    B_old = B;  
    V1DinvW = Nothing  
    #nzr21 = maximum([floor(0.6*nzr21),r])
    nzr_target = floor((1-w21)*nzr21 + w21*r)
    @show nzr_target
    while true
        V1DinvΦ = V1Dinv + Φ
        # update Z where W := V2*Z*U1'
        W  = V2V2T*(-V1DinvΦ + B)
        # update B
        F = V1DinvΦ + W
        svec2normW = sortperm(norm.(eachrow(F)),rev=true)
        B = zeros(n,r);
        for (i,idx) in enumerate(svec2normW)
            if i <= nzr_target 
                B[idx,:] = F[idx,:]
            else
                break
            end
        end
        # Compute residuals
        V1DinvW = V1Dinv + W
        res_infeasB = V1DinvW - B
        primal_resB = norm(res_infeasB)
        dual_resB = rho*norm(B-B_old)
        # Update Q
        Φ += res_infeasB
        # plots
        if iter >= 0
            append_to_plot!(P1,D1,It,primal_resB,dual_resB,iter)
        end
        # Stopping criteria
        iter += 1
        eps_pB = sqrt(n*r)*eps_abs + eps_rel*maximum([norm(B),norm(V1DinvW)])
        eps_dB = sqrt(n*r)*eps_abs + eps_rel*rho*norm(Φ)
        # if (primal_resP <= eps_pP && dual_resE <= eps_dE && primal_resQ <= eps_pQ && dual_resB <= eps_dB) || iter == max_iter
        if (primal_resB <= eps_pB && dual_resB <= eps_dB) || (iter >= max_iter)
            break
        end
        B_old = B
    end
    # Ψ = rho*Φ
    admmsol = SolutionADMM();
    admmsol.H = V1DinvW*U1';
    admmsol.iter = iter;
    admmsol.res_cp,admmsol.res_cd1,admmsol.res_cd2 = primal_resB,dual_resB,0.0;
    admmsol.eps_p,admmsol.eps_d = eps_pB,eps_dB/rho; 
    return admmsol
end

function getFig20norm(m1)
    global Pall,Dall,Itall = [],[],[];
    global D0all,D1all,D21all,D20all = [],[],[],[];
    global m = m1;
    global n,r = floor(Int64,0.5*m),floor(Int64,0.25*m);
    nameInst = string("A_",m,"_",n,"_",r);
    println(string("\nStarting instance: m,n,r: ",m,",",n,",",r))
    A = getMatlabInstance(nameInst,"A");
    inst = GinvInst(A,m,n,r);
    # Initialization
    U,S,V = svd(A,full=true)
    V1,V2 = V[:,1:r],V[:,r+1:n];
    U1,D = U[:,1:r],diagm(S[1:r]);
    
    V1Dinv,V1U1 = V1*inv(D),V1*U1';
    G = V1Dinv*U1';
    Λ1 = (1/maximum(abs.(V1U1)))*V1U1 
    Λ21 = (1/maximum(norm.(eachrow(V1))))*V1

    time_admm_21 = @elapsed admmsol_21 = runADMM21n(V1Dinv,V2,U1,Λ21,:ML,false);
    admmsol_21.z = getnorm21(admmsol_21.H);
    admmsol_21.time = time_admm_21;
    admmres_21 = getResultsADMM(inst,admmsol_21);
    #writeCSV!(arr_resADMM_21,admmres_21,Symbol(:For20_ADMM_21_,:OTM))
    println(string("m = ",m,". ADMM 2,1 finished in ", round_exact(admmsol_21.time,2), " sec. 2,1 norm ", round_exact(admmsol_21.z,3), ", 2,0 norm ", admmres_21.NZR, ". Iter: ", admmsol_21.iter));
    flush(stdout)

    for w_i = [25,50,75,80,90,95]
    # for w_i = [85,90,95,100]
        global P1,D0,D1,D21,D20,It = [],[],[],[],[],[];
        w21 = 0.01*w_i;
        time_admm_20 = @elapsed admmsol_20 = runADMM20n(V1Dinv,V2,U1,admmres_21.NZR,w21,true);
        push!(Pall,copy(P1));
        push!(D0all,copy(D0));
        push!(D1all,copy(D1));
        push!(D20all,copy(D20));
        push!(D21all,copy(D21));
        push!(Itall,copy(It));
        admmsol_20.z = getnorm21(admmsol_20.H);
        admmsol_20.time = time_admm_20;
        admmres_20 = getResultsADMM(inst,admmsol_20);
        #writeCSV!(arr_resADMM_20,admmres_20,Symbol(:ADMM_20_,w_i))
        println(string("w21 = ",w21,". ADMM 2,0 finished in ", round_exact(admmsol_20.time,2), " sec. 2,1 norm ", round_exact(admmsol_20.z,3), ", 2,0 norm ", admmres_20.NZR,". Iter: ", admmsol_20.iter));
        flush(stdout)
    end
    # save_plot(Itall,Pall,[".25" ".50" ".75" ".80" ".90" ".95"],"# iterations","primal residual gap",string("admm_20norm_primal6_m_",m1))
    save_plot20(Itall,D20all,[L"\omega = 0.25" L"\omega = 0.50" L"\omega = 0.75" L"\omega = 0.80" L"\omega = 0.90" L"\omega = 0.95"],"# iterations",L"{\|\|H\|\|}_{2,0}",string("admm_20norm_20n6_m_",m1))
    save_plot20(Itall,D0all,[L"\omega = 0.25" L"\omega = 0.50" L"\omega = 0.75" L"\omega = 0.80" L"\omega = 0.90" L"\omega = 0.95"],"# iterations",L"{\|\|H\|\|}_{0}",string("admm_20norm_0n6_m_",m1))
    save_plot20(Itall,D1all,[L"\omega = 0.25" L"\omega = 0.50" L"\omega = 0.75" L"\omega = 0.80" L"\omega = 0.90" L"\omega = 0.95"],"# iterations",L"{\|\|H\|\|}_{2,0}",string("26_1norm_20_m_",m1))
    save_plot20(Itall,D21all,[L"\omega = 0.25" L"\omega = 0.50" L"\omega = 0.75" L"\omega = 0.80" L"\omega = 0.90" L"\omega = 0.95"],"# iterations",L"{\|\|H\|\|}_{2,0}",string("26_21norm_20_m_",m1))
end


function getFigAllnorm(m1)
    global Pall,Dall,Itall = [],[],[];
    global Pall,Dall,Itall = [],[],[];
    global D0all,D1all,D21all,D20all = [],[],[],[];
    global m = m1;
    global n,r = floor(Int64,0.5*m),floor(Int64,0.25*m);
    nameInst = string("A_",m,"_",n,"_",r);
    println(string("\nStarting instance: m,n,r: ",m,",",n,",",r))
    A = getMatlabInstance(nameInst,"A");
    inst = GinvInst(A,m,n,r);
    # Initialization
    U,S,V = svd(A,full=true)
    V1,V2 = V[:,1:r],V[:,r+1:n];
    U1,D = U[:,1:r],diagm(S[1:r]);
    
    V1Dinv,V1U1 = V1*inv(D),V1*U1';
    G = V1Dinv*U1';
    Λ1 = (1/maximum(abs.(V1U1)))*V1U1 
    Λ21 = (1/maximum(norm.(eachrow(V1))))*V1

    global PR,DR,D0,D1,D20,D21,It = [],[],[],[],[],[],[];
    time_admm_1 = @elapsed admmsol_1 = runADMM1n(G,V2,U1,Λ1,:ML,true);
    admmsol_1.z = getnorm1(admmsol_1.H);
    admmsol_1.time = time_admm_1;
    admmres_1 = getResultsADMM(inst,admmsol_1);
    push!(Pall,copy(PR));
    push!(Dall,copy(DR));
    push!(D0all,copy(D0));
    push!(D1all,copy(D1));
    push!(D20all,copy(D20));
    push!(D21all,copy(D21));
    push!(Itall,copy(It));
    println(string("m = ",m,". ADMM 1 finished in ", round_exact(admmsol_1.time,2), " sec. 1 norm ", round_exact(admmsol_1.z,3), ", 2,0 norm ", admmres_1.NZR, ". Iter: ", admmsol_1.iter));
    flush(stdout)

    global PR,DR,D0,D1,D20,D21,It = [],[],[],[],[],[],[];
    time_admm_21 = @elapsed admmsol_21 = runADMM21n(V1Dinv,V2,U1,Λ21,:ML,true);
    admmsol_21.z = getnorm21(admmsol_21.H);
    admmsol_21.time = time_admm_21;
    admmres_21 = getResultsADMM(inst,admmsol_21);
    push!(Pall,copy(PR));
    push!(Dall,copy(DR));
    push!(D0all,copy(D0));
    push!(D1all,copy(D1));
    push!(D20all,copy(D20));
    push!(D21all,copy(D21));
    push!(Itall,copy(It));
    println(string("m = ",m,". ADMM 2,1 finished in ", round_exact(admmsol_21.time,2), " sec. 2,1 norm ", round_exact(admmsol_21.z,3), ", 2,0 norm ", admmres_21.NZR, ". Iter: ", admmsol_21.iter));
    flush(stdout)

    # global P1,D1,It = [],[],[];
    # time_admm_20 = @elapsed admmsol_20 = runADMM20nplot(V1Dinv,V2,U1,admmres_21.NZR,0.50);
    # push!(Pall,copy(P1));
    # push!(Dall,copy(D1));
    # push!(Itall,copy(It));
    # admmsol_20.z = getnorm21(admmsol_20.H);
    # admmsol_20.time = time_admm_20;
    # admmres_20 = getResultsADMM(inst,admmsol_20);
    # writeCSV!(arr_resADMM_20,admmres_20,Symbol(:ADMM_20_,0.50))
    # println(string("w21 = ",0.50,". ADMM 2,0 finished in ", round_exact(admmsol_20.time,2), " sec. 2,1 norm ", round_exact(admmsol_20.z,3), ", 2,0 norm ", admmres_20.NZR,". Iter: ", admmsol_20.iter));
    save_plot(Itall,Pall,[latexstring("ADMM\$_{1}\$") latexstring("ADMM\$_{2,1}\$")],"# iterations","primal residual",string("1vs21norm_pr_m_",m1))
    save_plot(Itall,Dall,[latexstring("ADMM\$_{1}\$") latexstring("ADMM\$_{2,1}\$")],"# iterations","dual residual",string("1vs21norm_dr_m_",m1))

    save_plot(Itall,D0all,[latexstring("ADMM\$_{1}\$") latexstring("ADMM\$_{2,1}\$")],"# iterations",L"{\|\|H\|\|}_{0}",string("1vs21norm_0n_m_",m1))
    save_plot(Itall,D1all,[latexstring("ADMM\$_{1}\$") latexstring("ADMM\$_{2,1}\$")],"# iterations",L"{\|\|H\|\|}_{1}",string("1vs21norm_1n_m_",m1))
    save_plot(Itall,D20all,[latexstring("ADMM\$_{1}\$") latexstring("ADMM\$_{2,1}\$")],"# iterations",L"{\|\|H\|\|}_{2,0}",string("1vs21norm_20n_m_",m1))
    save_plot(Itall,D21all,[latexstring("ADMM\$_{1}\$") latexstring("ADMM\$_{2,1}\$")],"# iterations",L"{\|\|H\|\|}_{2,1}",string("1vs21norm_21n_m_",m1))
end

global arr_resADMM_1,arr_resADMM_21,arr_resADMM_20,arr_resADMM_210 = init_ginv_res_admm(), init_ginv_res_admm(),init_ginv_res_admm(),init_ginv_res_admm();
# println("\nStarting procedure...")
M = [100,200,300,400,500,1000,2000,3000,4000,5000];
for m = [1000]
    getFigAllnorm(m)
end