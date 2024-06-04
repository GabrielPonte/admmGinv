function admm1norm(ginvInit::GinvInit;eps_abs=1e-4,eps_rel=1e-4,eps_opt=1e-5,rho=3,max_iter=1e5,time_limit=7200,stop_limit=:Boyd)
    # initialize timer
    time_start = time_ns()
    # initialize parameters
    U1,V2 = ginvInit.U1,ginvInit.V2;
    V1,V1DinvU1T = ginvInit.V1,ginvInit.V1DinvU1T;
    V2V2T = ginvInit.V2V2T;
    V1U1T,U1U1T = V1*U1',U1*U1';
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
    return admmsol,M,Λ
end

function admm2120norm(ginvInit::GinvInit,ω21::Float64,nzr21::Int64,M_21,Λ_21;eps_abs=1e-4,eps_rel=1e-4,rho=1e4,max_iter=1e5,time_limit=7200)
    # initialize timer
    time_start = time_ns()
    # initialize parameters
    U1,V2 = ginvInit.U1,ginvInit.V2;
    V1Dinv = ginvInit.V1Dinv;
    V2V2T = ginvInit.V2V2T;
    n,r = size(V1Dinv);
    Λ,E,M = Λ_21,M_21,M_21;
    E_old = copy(E)
    rho_inv = 1/rho; iter=0
    nzr_target = floor(Int64,(1-ω21)*nzr21 + ω21*r)
    c_val = 2; primal_res = Inf;
    # start admm
    while true
        rho_inv = 1/rho;
        # update Z where W := V2*J with J = (- V1DinvU1T + E - Λ)
        W = V2V2T*(E - rho_inv*Λ)
        M = V1Dinv + W
        # update of E
        E = closed_form2120n(M + rho_inv*Λ,rho_inv,nzr_target)
        res_infeas = M - E
        primal_res = norm(res_infeas)
        # Update P
        Λ += rho*res_infeas
        # Stopping criteria
        iter += 1
        n20 = getnorm20(M,1e-5)    
        if (iter == max_iter) || ((time_ns() - time_start)/1e9) >= 7200
            break
        end
        if n20 <= nzr_target && primal_res <= 1e-4
            dual_res = rho*norm(V2'*(E-E_old))
            eps_d = sqrt((n-r)*r)*eps_abs + eps_rel*norm(V2'*Λ)
            if dual_res <= eps_d
                break
            else 
                while c_val >= 1.2
                    W_try = V2V2T*(E - (c_val/rho)*Λ)
                    M_try = V1Dinv + W_try
                    E_try = closed_form2120n(M_try + (c_val/rho)*Λ,(c_val/rho),nzr_target)
                    if norm((M_try - E_try)) <= 1e-4
                        rho /= c_val
                        break
                    else
                        c_val -= 0.1
                    end
                end
            end
        end
        E_old = E;
    end
    admmsol = SolutionADMM();
    admmsol.time = (time_ns() - time_start)/1e9;
    admmsol.H = M*U1';
    admmsol.iter = iter;
    admmsol.res_pri,admmsol.res_dual,admmsol.res_d_ML = primal_res,rho*norm(V2'*(E-E_old)),rho*norm(V2'*(E-E_old));
    admmsol.z = getnorm21(admmsol.H);
    return admmsol
end