# function admm1norm(ginvInit::GinvInit;eps_abs=1e-4,eps_rel=1e-4,eps_opt=1e-4,rho=1e0,max_iter=1e3,time_limit=7200,stop_limit=:Boyd,adp=1)
#     # initialize timer
#     time_start = time_ns()
#     # initialize parameters
#     rho_inv = 1/rho; 
#     U1,V2 = ginvInit.U1,ginvInit.V2;
#     V1,V1DinvU1T = ginvInit.V1,ginvInit.V1DinvU1T;
#     U1U1T,V2V2T = U1*U1',ginvInit.V2V2T;

#     Θ0 = (1/maximum(abs.(V1*U1')))*V1*U1'
#     Θ1 = copy(Θ0); 

#     E0 = V1DinvU1T + Θ0; 
#     E1 = copy(E0); 


#     W0 = V2V2T*(E0 - rho_inv*Θ0)*U1U1T; 
#     W1 = copy(W0);
    
#     E = Nothing
#     H = Nothing; 
#     Θ = Nothing;
#     Θ_hat0 = Nothing;
#     Θ_hat1 = Nothing; 

#     n,m = size(Θ0);r = size(U1,2)
    
#     iter=1; norm_V1DinvU1T = norm(V1DinvU1T);
#     primal_res,dual_res,opt_res = 0.0,0.0,0.0;
#     eps_p,eps_d = 0.0,0.0;

#     pres = []; dres = []; mres = [];
#     rhos = [rho]; objs = []; tols = [];

#     #Λ = rho_inv*Θ;
#     # Fast ADMM, Nestorov
#     restart = 0.999;
#     alpha1 = 1;
#     # Residual balancing
#     bs = 2; rs = 0.1;
#     # Spectral
#     freq = 2;
#     minval = 1e-20;
#     orthval = max(0.2, minval);

#     tol = 1e-3; siter = 1; eiter = 1e5;

#     # start admm
#     while true
#         rho_inv = 1/rho;
#         # update Z where W := V2*J*U1' with J = (- V1DinvU1T + E - Λ)
#         W = V2V2T*(E1 - rho_inv*Θ1)*U1U1T
#         # update H
#         H = V1DinvU1T + W
#         # update of E
#         E = closed_form1n(H + rho_inv*Θ1,rho_inv)
#         # Compute residuals
#         pres1 = H - E;
#         dres1 = V2'*(E-E1)*U1;
#         # Update lagrangian variable
#         Θ = Θ1 + rho*pres1;

#         # Stopping criteria
#         push!(pres,norm(pres1)); push!(dres,rho*norm(dres1)); 
#         push!(mres,rho*pres[iter]^2+rho*norm(E-E1)^2);
#         push!(objs,getnorm1(H));

#         pres_norm = pres[iter]/max(norm(E),norm(W),norm_V1DinvU1T);
#         dres_norm = dres[iter]/norm(V2'*Θ*U1);
#         push!(tols,max(pres_norm, dres_norm));

#         if tols[iter] < tol || iter == max_iter
#             opt_res = abs(getnorm1(H)-tr(Θ'*V1DinvU1T))
#             break
#         end

#         if adp == 2 # Fast ADMM, Nesterov with restart

#             if iter == 1 || mres[iter] < restart*mres[iter-1]
#                 alpha = (1+sqrt(1+4*alpha1^2))/2;
#                 E1 = E + (alpha1-1)/alpha*(E-E0);
#                 Θ1 = Θ + (alpha1-1)/alpha*(Θ-Θ0);
#                 alpha1 = alpha;
#                 E0 = deepcopy(E);
#                 Θ0 = deepcopy(Θ);
#                 E1 = deepcopy(E); #previous Bv
#             else
#                 alpha = 1; alpha1 = 1;
#                 E1 = deepcopy(E0);
#                 Θ1 = deepcopy(Θ0);
#                 pres[iter] = pres[iter-1];
#                 dres[iter] = dres[iter-1];
#                 mres[iter] = mres[iter-1];
#             end

#         elseif adp == 3 # Residual balancing

#             if iter>siter && iter < eiter
#                 if dres[iter] < pres[iter] * rs #dual residual is smaller, need large rho
#                     rho = bs * rho;
#                 elseif pres[iter] < dres[iter] * rs #primal residual is smaller, need small rho
#                     rho = rho/bs;
#                 end
#             end 

#         elseif adp == 4 # Normalized residual balancing

#             if iter>siter && iter < eiter
#                 if dres_norm < pres_norm * rs #dual residual is smaller, need large rho
#                     rho = bs * rho;
#                 elseif pres_norm < dres_norm * rs #primal residual is smaller, need small rho
#                     rho = rho/bs;
#                 end
#             end 
            
#         elseif adp == 5  # Spectral penalty
#             if iter == 1
#                 Θ0 = deepcopy(Θ);
#                 Θ_hat0 = Θ1 + rho*(H-E1);
#                 E0 = deepcopy(E);
#                 W0 = deepcopy(W);
#             elseif mod(iter,freq)==0 && iter>siter && iter < eiter
#                 Θ_hat = Θ1 + rho*(H-E1);
#                 rho = adaptive_update(rho,W,W0,Θ_hat,Θ_hat0,E,E0,Θ,Θ0,orthval,minval)
#                 # record for next estimation
#                 Θ0 = deepcopy(Θ);
#                 Θ_hat0 = deepcopy(Θ_hat);
#                 E0 = deepcopy(E);
#                 W0 = deepcopy(W);
#             end
#         end
#         push!(rhos,rho);
#         E1 = deepcopy(E); Θ1 = deepcopy(Θ);iter += 1;

#     end
#     #plot(collect(1:length(tols)),tols)
#     #plot(collect(1:length(objs)),objs)
#     admmsol = SolutionADMM();
#     admmsol.time = (time_ns() - time_start)/1e9;
#     admmsol.H = H;
#     admmsol.iter = iter;
#     admmsol.res_pri,admmsol.res_dual,admmsol.res_d_ML = primal_res,norm(V2'*Θ*U1),rho*norm(V2'*(E-E1)*U1);
#     admmsol.res_opt = opt_res
#     admmsol.eps_p,admmsol.eps_d = eps_p,eps_d; 
#     admmsol.z = getnorm1(admmsol.H);
#     return admmsol,pres,dres,tols,objs,rhos
# end


function admm1norm(ginvInit::GinvInit;eps_abs=1e-4,eps_rel=1e-4,eps_opt=1e-5,rho=3,max_iter=1e5,time_limit=7200,stop_limit=:Boyd)
    # initialize timer
    time_start = time_ns()
    # initialize parameters
    U1,V2 = ginvInit.U1,ginvInit.V2;
    V1,V1DinvU1T = ginvInit.V1,ginvInit.V1DinvU1T;
    V2V2T = ginvInit.V2V2T;
    V1U1T,U1U1T = V1*U1',U1*U1';
    # Θ  = (1/maximum(abs.(V1U1T)))*V1U1T
    Θ = (1/norm(V1U1T,Inf))*V1U1T;
    m,r = size(U1)
    n  = size(V2,1)
    rho_inv = 1/rho;
    iter=0; norm_V1DinvU1T = norm(V1DinvU1T);
    primal_res,dual_res,opt_res = 0.0,0.0,0.0;
    eps_p,eps_d = 0.0,0.0;
    Λ = rho_inv*Θ; E = V1DinvU1T + Λ; 
    E_old = E; H = Nothing;
    while true
        # update Z where W := V2*Z*U1' with J = (- G + E - Λ)
        W = V2V2T*(E - Λ)*U1U1T
        # update E
        H = V1DinvU1T + W
        # update of E
        E = closed_form1n(H + Λ,rho_inv)
        # Compute residuals
        res_infeas = H - E
        primal_res = norm(res_infeas)
        dual_res = rho*norm(V2'*(E-E_old)*U1)
        # Update P
        Λ += res_infeas
        # Stopping criteria
        if stop_limit == :Boyd
            eps_p = sqrt(n*m)*eps_abs + eps_rel*maximum([norm(E),norm(W),norm_V1DinvU1T])
            eps_d = sqrt((n-r)*r)*eps_abs + eps_rel*rho*norm(V2'*Λ*U1)
        elseif stop_limit == :OptGap
            dual_res = rho*norm(V2'*Λ*U1)
            eps_p = eps_opt;eps_d = eps_opt
        else
            error("stop limit not defined correctly")
        end
        iter += 1
        if iter == max_iter || ((time_ns() - time_start)/1e9) >= 7200
            break
        end
        if (primal_res <= eps_p && dual_res <= eps_d)
            opt_res = abs(getnorm1(H)-rho*tr(V1DinvU1T'*Λ))
            break
        end
        E_old = E
    end
    
    admmsol = SolutionADMM();
    admmsol.time = (time_ns() - time_start)/1e9;
    E_diff = E-E_old
    admmsol.H = H;
    admmsol.iter = iter;
    admmsol.res_pri,admmsol.res_dual,admmsol.res_d_ML = primal_res,rho*norm(V2'*Λ*U1),rho*norm(V2'*E_diff*U1);
    admmsol.res_opt = opt_res
    admmsol.eps_p,admmsol.eps_d = eps_p,eps_d; 
    admmsol.z = getnorm1(admmsol.H);
    return admmsol
   
end

function admm21norm(ginvInit::GinvInit;eps_abs=1e-7,eps_rel=1e-7,eps_opt=1e-5,rho=1,max_iter=1e5,time_limit=7200,stop_limit=:Boyd)
    # initialize timer
    time_start = time_ns()
    # initialize parameters
    U1,V2 = ginvInit.U1,ginvInit.V2;
    V1,V1Dinv = ginvInit.V1,ginvInit.V1Dinv;
    V2V2T = ginvInit.V2V2T;
    Θ = (1/maximum(norm.(eachrow(V1))))*V1;
    n,r = size(Θ);
    rho_inv = 1/rho;
    iter=0; norm_V1Dinv = norm(V1Dinv);
    primal_res,dual_res,opt_res = 0.0,0.0,0.0;
    eps_p,eps_d = 0.0,0.0;
    Λ = rho_inv*Θ; E = ginvInit.V1Dinv+Λ; 
    E_old = E; M = Nothing;
    # start admm
    while true
        # update Z where W := V2*J with J = (- V1DinvU1T + E - Λ)
        W = V2V2T*(E - Λ)
        # update of E
        M = V1Dinv + W
        E = closed_form21n(M + Λ,rho_inv)
        # Compute residuals
        res_infeas = M - E
        # Update P
        Λ += res_infeas
        # Stopping criteria
        primal_res = norm(res_infeas)
        if stop_limit == :Boyd
            dual_res = rho*norm(V2'*(E-E_old))
            eps_p = sqrt(n*r)*eps_abs + eps_rel*maximum([norm(E),norm(W),norm_V1Dinv])
            eps_d = sqrt((n-r)*r)*eps_abs + eps_rel*rho*norm(V2'*Λ)
            E_old = E
        elseif stop_limit == :OptGap
            dual_res = rho*norm(V2'*Λ)
            eps_p = eps_opt,eps_d = eps_opt
        else
            error("stop limit not defined correctly")
        end
        iter += 1
        if (iter == max_iter) || ((time_ns() - time_start)/1e9) >= 7200 || (primal_res <= eps_p && dual_res <= eps_d)
            opt_res = abs(getnorm1(M)-rho*tr(Λ'*V1Dinv))
            break
        end
    end
    admmsol = SolutionADMM();
    admmsol.time = (time_ns() - time_start)/1e9;
    admmsol.H = M*U1';
    admmsol.iter = iter;
    admmsol.res_pri,admmsol.res_dual,admmsol.res_d_ML = primal_res,rho*norm(V2'*Λ),rho*norm(V2'*(E-E_old));
    admmsol.res_opt = opt_res
    admmsol.eps_p,admmsol.eps_d = eps_p,eps_d; 
    admmsol.z = getnorm21(admmsol.H);
    return admmsol
end

function admm20norm(ginvInit::GinvInit,ω21::Float64,nzr21::Int64;rho=1,max_iter=1e5,time_limit=7200)
    # initialize timer
    time_start = time_ns()
    # initialize parameters
    U1,V2 = ginvInit.U1,ginvInit.V2;
    V1Dinv = ginvInit.V1Dinv;
    V2V2T = ginvInit.V2V2T;
    #@show typeof(V1Dinv)
    n,r = size(V1Dinv);
    Φ,B = zeros(n,r),V1Dinv
    rho_inv = 1/rho; iter=0; M = Nothing;
    nzr_target = floor(Int64,(1-ω21)*nzr21 + ω21*r)
    # start admm
    while true
        # update Z where W := V2*J with J = (- V1DinvU1T + E - Λ)
        W = V2V2T*(B - Φ)
        # update of B
        M = V1Dinv + W
        B = closed_form20n(M + Φ,nzr_target)
        # Update P
        Φ += (M - B)
        # Stopping criteria
        iter += 1
        if (iter == max_iter) || ((time_ns() - time_start)/1e9) >= 7200 || (getnorm20(M,1e-5) <= nzr_target)
            break
        end
    end
    admmsol = SolutionADMM();
    admmsol.time = (time_ns() - time_start)/1e9;
    admmsol.H = M*U1';
    admmsol.iter = iter;
    admmsol.z = getnorm21(admmsol.H);
    return admmsol
end






function runADMM1n(G,V2,U1,Θ,TP::Symbol,PLT::Bool)
    # initial params
    time_start = time_ns()
    max_iter = 1e5;
    rho = 3;
    rho_inv = 1/rho
    V2V2T = V2*V2';
    U1U1T = U1*U1';
    iter = 0;
    primal_res,dual_res,opt_res = 0.0,0.0,0.0;
    Λ = rho_inv*Θ;
    E = G+Λ;
    E_old = E;    
    GW = Nothing
    norm_G = norm(G);
    if TP == :ML
        eps_abs,eps_rel = 1e-4,1e-4;
        eps_p,eps_d,eps_opt = 0.0,0.0,0.0
    else
        eps_p,eps_d,eps_opt = 1e-4,1e-4,1e-4
        
    end
    while true
        # update Z where W := V2*Z*U1' with J = (- G + E - Λ)
        W = V2V2T*(E - Λ)*U1U1T
        # update E
        GW = G + W
        Y = GW + Λ
        # update of E
        E = closed_form1n(Y,rho_inv)
        #E = sign.(Y).*max.(abs.(Y) .- (rho_inv), 0)
        # E2 = subprobE(n,Y,rho)
        # @show norm(E2-E)
        # flush(stdout)
        # Compute residuals
        E_diff = E-E_old
        res_infeas = GW - E
        primal_res = norm(res_infeas)
        dual_res = rho*norm(V2'*E_diff*U1)
        # Update P
        Λ += res_infeas
        # Stopping criteria
        if TP == :ML
            eps_p = sqrt((n-r)*r)*eps_abs + eps_rel*maximum([norm(E),norm(W),norm_G])
            eps_d = sqrt((n-r)*r)*eps_abs + eps_rel*rho*norm(V2'*Λ*U1)
        end
        # else
        #     dual_res = rho*norm(V2'*Λ*U1)
        # end
        # plots
        # PLT ? append_to_plot!(P1,D1,It,primal_res,dual_res,iter) : "nvm"
        if PLT
            H = GW
            push!(PR,norm(res_infeas))
            push!(DR,norm(dual_res))
            push!(D20,getnorm20(H,1e-5))
            push!(D21,getnorm21(H))
            push!(D0,getnorm0(H,1e-5))
            push!(D1,getnorm1(H))
            push!(It,iter)
        end
        iter += 1
        if iter == max_iter || ((time_ns() - time_start)/1e9) >= 7200
            break
        end
        if (primal_res <= eps_p && dual_res <= eps_d)
            opt_res = abs(getnorm1(GW)-rho*tr(G'*Λ))
            break
            # if TP == :ML || opt_res <= eps_opt
            #     break
            # end
        end
        E_old = E
    end
    E_diff = E-E_old
    admmsol = SolutionADMM();
    admmsol.H = GW;
    admmsol.iter = iter;
    admmsol.res_pri,admmsol.res_dual,admmsol.res_d_ML = primal_res,rho*norm(V2'*Λ*U1),rho*norm(V2'*E_diff*U1);
    admmsol.res_opt = opt_res
    admmsol.eps_p,admmsol.eps_d = eps_p,eps_d; 
    return admmsol
end

function runADMM21n(V1Dinv,V2,U1,Θ,TP::Symbol,PLT::Bool)
    # initial params
    max_iter = 5e4;
    rho = 1;
    rho_inv = 1/rho
    iter = 0;
    primal_res,dual_res,opt_res = 0.0,0.0,0.0;
    V2V2T = V2*V2'
    Λ = rho_inv*Θ;
    E = V1Dinv + Λ;
    E_old = E;  
    GW = Nothing  
    normV1Dinv = norm(V1Dinv);
    if TP == :ML
        eps_abs,eps_rel = 1e-7,1e-7;
        eps_p,eps_d,eps_opt = 0.0,0.0,0.0
    else
        eps_p,eps_d,eps_opt = 1e-5,1e-5,1e-5
        
    end
    while true
        V1DinvΛ = V1Dinv + Λ
        # update Z where W := V2*Z
        W  = V2V2T*(E-Λ)
        # update E
        Y = V1DinvΛ + W
        E = closed_form21n(Y,rho_inv)
        
        for (i,el) in enumerate(norm_rows_y)
            if rho_inv < el
                E[i,:] = ((el - rho_inv)/el) * Y[i,:]
            end
        end
        # Compute residuals
        E_diff = E-E_old
        GW = V1Dinv + W
        res_infeas = GW -E
        primal_res = norm(res_infeas)
        dual_res = rho*norm(V2'*E_diff)
        # Update P
        Λ += res_infeas
        # Stopping criteria
        if TP == :ML
            eps_p = sqrt((n-r)*r)*eps_abs + eps_rel*maximum([norm(E),norm(W),normV1Dinv])
            eps_d = sqrt((n-r)*r)*eps_abs + eps_rel*rho*norm(V2'*Λ)
        end
        # else
        #     dual_res = rho*norm(V2'*Λ)
        # end
        # plots
        # PLT ? append_to_plot!(P1,D1,It,primal_res,dual_res,iter) : "nvm"
        if PLT
            H = GW*U1'
            push!(PR,norm(res_infeas))
            push!(DR,norm(dual_res))
            push!(D20,getnorm20(H,1e-5))
            push!(D21,getnorm21(H))
            push!(D0,getnorm0(H,1e-5))
            push!(D1,getnorm1(H))
            push!(It,iter)
        end
        iter += 1
        if iter == max_iter
            break
        end
        if (primal_res <= eps_p && dual_res <= eps_d)
            opt_res = abs(getnorm21(GW)-rho*tr(V1Dinv'*Λ))
            break
            # if TP == :ML || opt_res <= eps_opt
            #     break
            # end
        end
        E_old = E
    end
    E_diff = E-E_old
    admmsol = SolutionADMM();
    admmsol.H = GW*U1';
    admmsol.iter = iter;
    admmsol.res_pri,admmsol.res_dual,admmsol.res_d_ML = primal_res,rho*norm(V2'*Λ),rho*norm(V2'*E_diff);
    admmsol.res_opt = opt_res
    admmsol.eps_p,admmsol.eps_d = eps_p,eps_d; 
    return admmsol
end

function runADMM20n(V1Dinv,V2,U1,nzr21,w21::Float64,PLT::Bool)
    # initial params
    time_start = time_ns()
    max_iter = 5e4;
    rho = 1;
    rho_inv = 1/rho
    # eps_abs,eps_rel = 1e-7,1e-7*rho_inv;
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
        # 
        # 
        # Update Q
        Φ += res_infeasB
        # plots
        # append_to_plot!(P1,D1,It,primal_resQ,dual_resB,iter)
        # Stopping criteria
        norm20H = getnorm20(V1DinvW,1e-5)
        # plots
        if PLT
            H = V1DinvW*U1'
            push!(P1,norm(res_infeasB))
            push!(D20,getnorm20(H,1e-5))
            push!(D21,getnorm21(H))
            push!(D0,getnorm0(H,1e-5))
            push!(D1,getnorm1(H))
            push!(It,iter)
        end
        iter += 1
        if iter == max_iter || ((time_ns() - time_start)/1e9) >= 7200
            primal_resB = norm(res_infeasB)
            break
        end
        if norm20H <= nzr_target
            primal_resB = norm(res_infeasB)
            break
        end
        # eps_pB = sqrt(n*r)*eps_abs + eps_rel*maximum([norm(B),norm(V1DinvW)])
        # eps_dB = sqrt(n*r)*eps_abs + eps_rel*rho*norm(Φ)
        # # if (primal_resP <= eps_pP && dual_resE <= eps_dE && primal_resQ <= eps_pQ && dual_resB <= eps_dB) || iter == max_iter
        # if (primal_resB <= eps_pB && dual_resB <= eps_dB) || (iter >= max_iter)
        #     break
        # end
        B_old = B
    end
    # Ψ = rho*Φ
    admmsol = SolutionADMM();
    admmsol.H = V1DinvW*U1';
    admmsol.iter = iter;
    B_diff = B -B_old
    admmsol.res_pri,admmsol.res_dual,admmsol.res_d_ML = primal_resB,rho*norm(V2'*Φ),rho*norm(V2'*B_diff);
    admmsol.res_opt = 0.0
    admmsol.eps_p,admmsol.eps_d = 1e-5,1e-5; 
    return admmsol
end

function runADMM210nv3(V1Dinv,V2,U1,Θ,nzr_target)
    TP = :ML
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
        eps_abs,eps_rel = 1e-9,1e-9;
        eps_p,eps_d = 0.0,0.0
    else
        eps_p,eps_d = 1e-7,1e-7*rho_inv
    end
    @show nzr_target
    while true
       #cs = iter/max_iter
        #i_max = floor(Int64,cs*nzr_target + (1-cs*n))
        V1DinvΛ = V1Dinv + Λ
        # update Z where W := V2*Z
        W  = V2V2T*(-V1DinvΛ + E)
        # update E
        Y = V1DinvΛ + W
        E = zeros(n,r)
        norm_rows_y = norm.(eachrow(Y))
        svec2normW = sortperm(norm_rows_y,rev=true)
        for (i,idx) in enumerate(svec2normW)
            # if i <= nzr_target
            el = norm_rows_y[idx]
            if rho_inv < el
                E[i,:] = ((el - rho_inv)/el) * Y[i,:]
            end
            # else
            #     break
            # end
        end
        # Compute residuals
        GW = V1Dinv + W
        res_infeas = GW -E
        primal_res = norm(res_infeas)
        # Update P
        Λ += res_infeas
        # plots
        # append_to_plot!(P1,D1,It,primal_res,dual_res,iter)
        # Stopping criteria
        iter += 1
        if TP == :ML
            eps_p = sqrt(n*r)*eps_abs + eps_rel*maximum([norm(E),norm(GW)])
            eps_d = sqrt(n*r)*eps_abs + eps_rel*rho*norm(Λ)
            dual_res = rho*norm(E-E_old)
        else
            dual_res = norm(V2'*Λ)
        end
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

function runADMM210nv2(V1Dinv,V2,U1,Θ,nzr21)
    # initial params
    max_iter = 1000;
    eps_abs,eps_rel = 1e-4,1e-4;
    rho = 1e-6;
    rho_inv = 1/rho
    V2V2T = V2*V2'
    Λ = rho_inv*Θ
    E = V1Dinv + Λ
    Φ,B = zeros(n,r),V1Dinv
    W = zeros(n,r)
    iter = 0;
    primal_resE,dual_resE = 0.0,0.0;
    primal_resB,dual_resB = 0.0,0.0;
    eps_pB,eps_dB = 0.0,0.0;
    eps_pE,eps_dE = 0.0,0.0;
    E_old,B_old = E,B;  
    V1DinvW = Nothing  
    # w21 = 0.7;
    # nzr_target = floor((w21*nzr21 + (1-w21)*r))
    nzr_target = nzr21
    @show nzr_target
    while true
        V1DinvΦ = V1Dinv + Φ
        # V1DinvΛ = V1Dinv + Λ


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
        # update E
        # Y = V1DinvΛ + W
        Y = B - Φ
        E = zeros(n,r)
        norm_rows_y = norm.(eachrow(Y))
        for (i,el) in enumerate(norm_rows_y)
            if rho_inv < el
                E[i,:] = ((el - rho_inv)/el) * Y[i,:]
            end
        end
        # update Z where W := V2*Z*U1'
        # W  = 0.5*V2V2T*(B + E - V1DinvΦ - V1DinvΛ)
        W = E - V1Dinv
        # Compute residuals
        V1DinvW = V1Dinv + W
        res_infeasB = V1DinvW - B
        # res_infeasE = V1DinvW - E
        primal_res = norm([res_infeasB])
        # primal_resE = norm(res_infeasE)
        dual_res = rho*norm([B-B_old])
        # @show primal_res,dual_res
        # dual_resE = rho*norm(E-E_old)
        # Update Q
        Φ += res_infeasB
        # Λ += res_infeasE
        # plots
        # append_to_plot!(P1,D1,It,primal_resQ,dual_resB,iter)
        append_to_plot_extended!(primal_res,primal_res,dual_res,dual_res,iter)
        # Stopping criteria
        iter += 1
        eps_pB = sqrt(n*r)*eps_abs + eps_rel*maximum([norm([B]),norm([V1DinvW])])
        eps_dB = sqrt(n*r)*eps_abs + eps_rel*rho*norm([Φ])
        #eps_pB,eps_dB = Inf,Inf
        # eps_pE = sqrt(n*r)*eps_abs + eps_rel*maximum([norm(E),norm(V1DinvW)])
        # eps_dE = sqrt(n*r)*eps_abs + eps_rel*rho*norm(Λ)
        # if (primal_resB <= eps_pB && dual_resB <= eps_dB && primal_resE <= eps_pE && dual_resE <= eps_dE) || iter == max_iter
        if (primal_res <= eps_pB && dual_res <= eps_dB) || (iter >= max_iter)
            @show primal_res,dual_res,eps_pB,eps_dB,iter
            break
        end
        E_old = copy(E)
        # B_old = B
        # update rho
        # if primal_res >= 1.50*dual_res
        #     rho *= 2
        #     Φ,Λ = 0.5*Φ,0.5*Λ
        # elseif dual_res >= 1.50*primal_res
        #     rho *= 0.5
        #     Φ,Λ = 2*Φ,2*Λ
        # end
    end
    # Ψ = rho*Φ
    admmsol = SolutionExtendedADMM(); 
    admmsol.H = V1DinvW*U1';
    admmsol.iter = iter;
    admmsol.res_cpB,admmsol.res_cdB = primal_resB,dual_resB;
    admmsol.res_cpE,admmsol.res_cdE = primal_resE,dual_resE;
    admmsol.eps_pB,admmsol.eps_dB = eps_pB,eps_dB/rho; 
    admmsol.eps_pE,admmsol.eps_dE = eps_pE,eps_dE/rho; 
    return admmsol
end

function runADMM210n(V1Dinv,V2,U1,Θ,nzr21)
    # initial params
    max_iter = 5000;
    eps_abs,eps_rel = 1e-6,1e-6;
    rho = 1;
    rho_inv = 1/rho
    V2V2T = V2*V2'
    Λ = rho_inv*Θ
    E = V1Dinv + Λ
    Φ,B = zeros(n,r),V1Dinv
    iter = 0;
    primal_resE,dual_resE = 0.0,0.0;
    primal_resB,dual_resB = 0.0,0.0;
    eps_pB,eps_dB = 0.0,0.0;
    eps_pE,eps_dE = 0.0,0.0;
    E_old,B_old = E,B;  
    V1DinvW = Nothing  
    # w21 = 0.7;
    nzr_target = nzr21
    @show nzr_target
    while true
        V1DinvΦ = V1Dinv + Φ
        V1DinvΛ = V1Dinv + Λ
        # update Z where W := V2*Z*U1'
        W  = 0.5*V2V2T*(B + E - V1DinvΦ - V1DinvΛ)
        #W  = V2V2T*( E  - V1DinvΛ)
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
        V1DinvW = V1Dinv + W
        res_infeasB = V1DinvW - B
        res_infeasE = V1DinvW - E
        primal_res = norm([res_infeasB;res_infeasE])
        # primal_resE = norm(res_infeasE)
        dual_res = rho*norm([B-B_old;E-E_old])
        # @show primal_res,dual_res
        # dual_resE = rho*norm(E-E_old)
        # Update Q
        Φ += res_infeasB
        Λ += res_infeasE
        # plots
        # append_to_plot!(P1,D1,It,primal_resQ,dual_resB,iter)
        append_to_plot_extended!(primal_res,primal_res,dual_res,dual_res,iter)
        # Stopping criteria
        iter += 1
        eps_pB = sqrt(2*n*r)*eps_abs + eps_rel*maximum([norm([B;E]),norm([V1DinvW;V1DinvW])])
        eps_dB = sqrt(2*n*r)*eps_abs + eps_rel*rho*norm([Φ;Λ])
        #eps_pB,eps_dB = Inf,Inf
        # eps_pE = sqrt(n*r)*eps_abs + eps_rel*maximum([norm(E),norm(V1DinvW)])
        # eps_dE = sqrt(n*r)*eps_abs + eps_rel*rho*norm(Λ)
        # if (primal_resB <= eps_pB && dual_resB <= eps_dB && primal_resE <= eps_pE && dual_resE <= eps_dE) || iter == max_iter
        if (primal_res <= eps_pB && dual_res <= eps_dB) || (iter >= max_iter)
            @show primal_res,dual_res,eps_pB,eps_dB,iter
            break
        end
        E_old = E
        B_old = B
        # update rho
        # if primal_res >= 1.50*dual_res
        #     rho *= 2
        #     Φ,Λ = 0.5*Φ,0.5*Λ
        # elseif dual_res >= 1.50*primal_res
        #     rho *= 0.5
        #     Φ,Λ = 2*Φ,2*Λ
        # end
    end
    # Ψ = rho*Φ
    admmsol = SolutionExtendedADMM(); 
    admmsol.H = V1DinvW*U1';
    admmsol.iter = iter;
    admmsol.res_cpB,admmsol.res_cdB = primal_resB,dual_resB;
    admmsol.res_cpE,admmsol.res_cdE = primal_resE,dual_resE;
    admmsol.eps_pB,admmsol.eps_dB = eps_pB,eps_dB/rho; 
    admmsol.eps_pE,admmsol.eps_dE = eps_pE,eps_dE/rho; 
    return admmsol
end


# function runADMM20nFull(G,V2,U1,Φ)
#     # initial params
#     max_iter = 2500;
#     eps_abs,eps_rel = 1e-5,1e-5;
#     rho = 1;
#     rho_inv = 1/rho
#     iter = 0;
#     V2V2T = V2*V2';
#     U1U1T = U1*U1';
#     primal_resQ,dual_resB = 0.0,0.0;
#     Q,B = rho_inv*Φ,G;
#     B_old = B;  
#     GW = Nothing  
#     while true
#         GQ = G + Q
#         # update Z where W := V2*Z*U1'
#         W  = V2V2T*(-GQ + B)*U1U1T
#         # update B
#         F = GQ + W
#         svec2normW = sortperm(norm.(eachrow(F)),rev=true)
#         B = zeros(n,m);
#         for (i,idx) in enumerate(svec2normW)
#             if i <= r 
#                 B[idx,:] = F[idx,:]
#             else
#                 break
#             end
#         end
#         # Compute residuals
#         GW = G + W
#         res_infeasQ = GW - B
#         primal_resQ = norm(res_infeasQ)
#         dual_resB = rho*norm(B-B_old)
#         # Update Q
#         Q += res_infeasQ
#         # plots
#         append_to_plot!(P1,D1,It,primal_resQ,dual_resB,iter)
#         # Stopping criteria
#         iter += 1
#         eps_pQ = sqrt(n*m)*eps_abs + eps_rel*maximum([norm(B),norm(GW)])
#         eps_dB = sqrt(n*m)*eps_abs + eps_rel*rho*norm(Q)
#         # if (primal_resP <= eps_pP && dual_resE <= eps_dE && primal_resQ <= eps_pQ && dual_resB <= eps_dB) || iter == max_iter
#         if (primal_resQ <= eps_pQ && dual_resB <= eps_dB) || iter == max_iter

#             # @show primal_res,dual_res,iter
#             break
#         end
#         B_old = B
#     end
#     Φ = rho*Q
#     admmsol = SolutionADMM();
#     admmsol.H = GW;
#     admmsol.iter = iter;
#     admmsol.iter_p,admmsol.iter_d,admmsol.iter_u = 0,0,iter;
#     admmsol.res_cp,admmsol.res_cd1,admmsol.res_cd2 = primal_resQ,dual_resB,0.0;
#     return admmsol
# end

# function runADMM210n(G1,V2,U1)
#     # initial params
#     Z = Nothing
#     Λ,Φ = zeros(n,r), zeros(n,r)
#     max_iter = 1500;
#     eps_abs,eps_rel = 1e-5,1e-4;
#     rho = 1;
#     rho_inv = 1/rho
#     iter = 0;
#     primal_res,dual_res = 0.0,0.0;
#     P = rho_inv*Λ;
#     Q = rho_inv*Φ;
#     E,B = G1,G1;
#     E_old,B_old = E,B;    
#     while true
#         # V1DinvP = V1Dinv + P
#         GQ = G1 + Q
#         # update Z 
#         # Z = 0.5*V2'*(-V1DinvP - V1DinvQ + E + B)
#         Z = V2'*( -GQ + B)
#         # update E
#         V2Z = V2*Z;
#         # Y = V1DinvP + V2Z
#         # E = zeros(n,r)
#         # norm_rows_y = norm.(eachrow(Y))
#         # for (i,el) in enumerate(norm_rows_y)
#         #     if rho_inv < el
#         #         E[i,:] = ((el - rho_inv)/el) * Y[i,:]
#         #     end
#         # end
#         # update B
#         W = GQ + V2Z;
#         svec2normW = sortperm(norm.(eachrow(W)),rev=true)
#         B = zeros(n,r);
#         for (i,idx) in enumerate(svec2normW)
#             if i <= r 
#                 B[idx,:] = W[idx,:]
#             else
#                 break
#             end
#         end
#         # Compute residuals
#         GW = G1 + V2Z
#         # res_infeasP = GW - E
#         res_infeasQ = GW - B
#         # primal_resP = norm(res_infeasP)
#         primal_resQ = norm(res_infeasQ)
#         # dual_resE = rho*norm(E-E_old)
#         dual_resB = rho*norm(B-B_old)
#         # Update P
#         # P += res_infeasP
#         Q += res_infeasQ
#         # plots
#         append_to_plot!(P1,D1,It,primal_resQ,dual_resB,iter)
#         # append_to_plot_extended!(P1,D1,It,primal_resP,primal_resQ,dual_resE,dual_resB,iter)
#         # Stopping criteria
#         iter += 1
#         # eps_pP = sqrt(n*r)*eps_abs + eps_rel*maximum([norm(E),norm(GW)])
#         # eps_dE = sqrt(n*r)*eps_abs + eps_rel*rho*norm(P)
#         eps_pQ = sqrt(n*m)*eps_abs + eps_rel*maximum([norm(B),norm(GW)])
#         eps_dB = sqrt(n*m)*eps_abs + eps_rel*rho*norm(Q)
#         # if (primal_resP <= eps_pP && dual_resE <= eps_dE && primal_resQ <= eps_pQ && dual_resB <= eps_dB) || iter == max_iter
#         if (primal_resQ <= eps_pQ && dual_resB <= eps_dB) || iter == max_iter

#             # @show primal_res,dual_res,iter
#             break
#         end
#         # E_old = E
#         B_old = B
#     end
#     Λ = rho*Q
#     admmsol = SolutionADMM();
#     @show size(G1 + V2*Z)
#     admmsol.H = (G1 + V2*Z)*U1';
#     admmsol.iter = iter;
#     admmsol.iter_p,admmsol.iter_d,admmsol.iter_u = 0,0,iter;
#     admmsol.res_cp,admmsol.res_cd1,admmsol.res_cd2 = primal_res,norm(V2'*Λ),maximum(norm.(eachrow(Λ)))-1;
#     return admmsol
# end