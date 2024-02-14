function update_E21n(V1DinvΛ,V2Z,mu)
    Y = V1DinvΛ + V2Z
    E = zeros(n,r)
    norm_rows_y = norm.(eachrow(Y))
    for (i,el) in enumerate(norm_rows_y)
        if (1/mu) < el
            E[i,:] = ((el - (1/mu))/el) * Y[i,:]
        end
    end
    return E
end

function update_Z21n(V1DinvΛ,V2,E)
    return V2'*(-V1DinvΛ + E)
end

function runADMM1n(G,V2,U1,Λ)
    max_iter = 500;
    global mu = 2
    global Y = zeros(n,m)
    W = zeros(n,m)
    V2V2T = V2*V2';
    U1U1T = U1*U1';
    break_error = false
    iter = 0;
    rho_p,rho_d,rho_u = 1.1,1.005,1.025;
    iter_p,iter_d,iter_u = 0,0,0;
    primal_res,dual_res = 0.0,0.0;
    while true
        # Update E
        mu_inv = 1/mu
        Y = G + W + mu_inv*Λ
        E = sign.(Y).*max.(abs.(Y) .- (mu_inv), 0)
        # update Z where W := V2*Z*U1'
        J = (- G + E - (1/mu)*Λ);
        W = V2V2T*J*U1U1T 
        # Compute residuals
        res_infeas = G + W - E
        primal_res = norm(res_infeas)
        dual_res = norm(Λ,Inf)-1 #maximum(abs.(Λ))-1
        # Update Λ
        Λ += mu*res_infeas
        mu = 1.04*mu
        # Update mu
        # if primal_res >10*dual_res
        #     mu = rho_p*mu
        #     iter_p += 1
        # elseif dual_res > 10*primal_res
        #     mu = rho_d*mu
        #     iter_d += 1
        # else
        #     mu = rho_u*mu
        #     iter_u += 1
        # end
        
        iter += 1
        if iter%100 == 0
            @show primal_res,dual_res,iter,mu
           # @show 
        end
        if  (primal_res < 1e-4 && dual_res <= 0.1)
            rho_d = 1.01
        end
        if (primal_res < 1e-4 && dual_res <= 0.01) || iter == max_iter
            if !(norm(V2'*Λ*U1) >= 1e-5 && iter < max_iter)
                @show primal_res,dual_res,iter,mu
                break
            else
                break_error = true
            end
        end
    end
    # mu = 30
    


    admmsol = SolutionADMM();
    admmsol.H = G + W;
    admmsol.iter = iter;
    admmsol.iter_p,admmsol.iter_d,admmsol.iter_u = iter_p,iter_d,iter_u;
    admmsol.res_cp,admmsol.res_cd1,admmsol.res_cd2 = primal_res,norm(V2'*Λ*U1),dual_res;
    admmsol.break_error = break_error
    return admmsol
end



function runADMM21n(V1Dinv,V2,U1,Λ)
    max_iter = 5000
    mu = 1e-3
    Z = Nothing
    V2Z = zeros(n,r);
    break_error = false
    iter = 0;
    rho_p,rho_d,rho_u = 1.1,1.005,1.01;
    iter_p,iter_d,iter_u = 0,0,0;
    primal_res,dual_res = 0.0,0.0;
    while true
        V1DinvΛ = V1Dinv + (1/mu)*Λ
        E = update_E21n(V1DinvΛ,V2Z,mu)
        Z = update_Z21n(V1DinvΛ,V2,E)
        V2Z = V2*Z;
        res_infeas = V1Dinv + V2Z-E
        # Compute residuals
        primal_res = norm(res_infeas)
        dual_res = maximum(norm.(eachrow(Λ)))-1
        # Update Λ
        Λ += mu*res_infeas
        # Update mu
        if primal_res > 10*dual_res
            mu = rho_p*mu 
            iter_p += 1
        elseif dual_res > 10*primal_res
            mu = rho_d*mu
            iter_d += 1
        else
            mu = rho_u*mu
            iter_u += 1
        end
        iter += 1
        if (primal_res < 1e-5 && dual_res <= 0.01) || iter == max_iter
            if !(norm(V2'*Λ) >= 1e-5 && iter < max_iter)
                break
            else
                break_error = true
            end
        end
    end
    admmsol = SolutionADMM();
    admmsol.H = V1Dinv*U1' +V2*Z*U1';
    admmsol.iter = iter;
    admmsol.iter_p,admmsol.iter_d,admmsol.iter_u = iter_p,iter_d,iter_u;
    admmsol.res_cp,admmsol.res_cd1,admmsol.res_cd2 = primal_res,norm(V2'*Λ),dual_res;
    admmsol.break_error = break_error
    return admmsol
end