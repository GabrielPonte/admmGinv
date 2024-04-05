function curv_adaptive_BB(al_h, de_h)
    # adapive BB, reference: FASTA paper
    tmph = de_h/al_h; # correlation
    if tmph > .5
        tau_h = de_h;
    else
        tau_h = al_h - 0.5*de_h;
    end
    return tau_h
end

function adaptive_update(rho,W,W0,Θ_hat,Θ_hat0,E,E0,Θ,Θ0,orthval,minval)
    
    # inner product
    tmp = real(conj(W-W0).*(Θ_hat-Θ_hat0));
    ul_hat = sum(tmp);
    tmp = real(conj(E-E0).*(Θ-Θ0));
    vl = sum(tmp);

    # norm of lambda, lambda_hat
    tmp = Θ_hat-Θ_hat0;
    dl_hat = norm(tmp);
    tmp = Θ-Θ0;
    dl = norm(tmp);
    
    # norm of gradient change
    tmp = W-W0;
    du = norm(tmp);
    tmp = E-E0;
    dv = norm(tmp);
    # flag to indicate whether the curvature can be estimated
    hflag = false;
    gflag = false;
    # estimate curvature, only if it can be estimated
    # use correlation/othogonal to test whether can be estimated
    if ul_hat > orthval*du*dl_hat + minval
        hflag = true;
        al_h = dl_hat^2/ul_hat;
        de_h = ul_hat/du^2;
        bb_h = curv_adaptive_BB(al_h, de_h);
    end
    if vl > orthval*dv*dl + minval
        gflag = true;
        al_g = dl^2/vl;
        de_g = vl/dv^2;
        bb_g = curv_adaptive_BB(al_g, de_g);
    end
    #if curvature can be estimated for both terms, balance the two
    #if one of the curvature cannot be estimated, use the estimated one
    #or use the previous stepsize to estimate
    if hflag && gflag
        rho = sqrt(bb_h*bb_g);
    elseif hflag
        rho = bb_h;
    elseif gflag
        rho = bb_g;
    end
    return rho



    # Θ_hat0 = copy(Θ_hat);
    # Θ_hat  = Θ + rho*(M-E_old);
    # # get Δ
    # ΔΘ = (Θ - Θ0);
    # ΔΘ_hat = (Θ_hat - Θ_hat0); 
    # ΔW = (W_old - W0);
    # ΔE = -(E_old - E0);
    # δθ = norm(ΔΘ)^2; 
    # δθh = norm(ΔΘ_hat)^2; 
    # δw = norm(ΔW)^2;
    # δe = norm(ΔE)^2;
    # δθw = tr(ΔΘ'*ΔW);#(dot(vec(ΔΘ),vec(ΔW)));
    # δθe = tr(ΔΘ'*ΔE);#(dot(vec(ΔΘ),vec(ΔE))); #tr(ΔΘ'*ΔE);
    # δθhw = tr(ΔΘ_hat'*ΔW);#(dot(vec(ΔΘ),vec(ΔW)));
    # δθhe = tr(ΔΘ_hat'*ΔE);#(dot(vec(ΔΘ),vec(ΔE))); #tr(ΔΘ'*ΔE);
    # # estimate spectral alphas
    # alpha_sd = δθh/δθhw;
    # alpha_mg = δθhw/δw;
    # if 2*alpha_mg > alpha_sd
    #     alpha = alpha_mg;
    # else
    #     alpha = alpha_sd - alpha_mg/2;
    # end
    # # estimate spectral betas
    # beta_sd = δθ/δθe;
    # beta_mg = δθe/δe;
    # if 2*beta_mg > beta_sd
    #     beta = beta_mg;
    # else
    #     beta = beta_sd - beta_mg/2;
    # end
    
    # # estimate correlations
    # alpha_cor = δθw/(norm(ΔΘ)*norm(ΔW));
    # beta_cor  = δθe/(norm(ΔΘ)*norm(ΔE));
    # # update rho
    # ϵcor = 0.2;
    # if alpha_cor > ϵcor && beta_cor > ϵcor
    #     rho_new = sqrt(alpha*beta);
    # elseif alpha_cor > ϵcor && beta_cor <= ϵcor
    #     rho_new = alpha;
    # elseif alpha_cor <= ϵcor && beta_cor > ϵcor
    #     rho_new = beta;
    # else 
    #     rho_new = rho;
    # end
    # # if isnan(alpha) || isnan(beta)
    # #     @show δθw,δθe,δθ,δw,δθw
    # #     @show alpha,beta
    # #     @show rho,rho_new
    # #     error("..")
    # # end
    # # @show alpha_cor,beta_cor,ϵcor
    # return rho_new,Θ_new
end