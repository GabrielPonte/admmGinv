getLICols(A,r) = (qr(A, ColumnNorm()).p)[1:r]
getLIRows(A,r) = (qr(A', ColumnNorm()).p)[1:r]

function get_LI_RC(A,r)
    time_start = time_ns() 
    C = getLICols(A,r);
    R = getLIRows(A[:,C],r);
    min_svd = minimum(svd(A[R,C]).S);
    if min_svd < 1e-5
        @warn "Submatrix close to be singular. Minimum singular value is $(@sprintf("%.2e", min_svd))";
    end
    elapsed_time = (time_ns() - time_start)/1e9
    return R,C,elapsed_time
end
    
function lu_rank_one_update(L, U, y, z, r)
    # only works for P = (1:r) otherwise y = y[P]
    # z = ei
    # y = A(R,C(i)) - A(R,C(alpha_idx))
    # y = y[P]; 
    for i = (1:r)
        U[i, i] = U[i, i] + y[i]*z[i];
        z[i] = z[i] / U[i, i]; 
        for j = (i+1:r)
            y[j] = y[j] - y[i]*L[j,i];
            L[j,i] = L[j,i] + z[i]*y[j];
            U[i,j] = U[i,j] + y[i]*z[j];
            z[j] = z[j] - z[i]*U[i,j];
        end
    end
    return L,U
end


function LS_det(A,R,C,TP::Symbol)
    time_start = time_ns()
    m,n = size(A); r = length(C);
    eps = 1e-4;
    swaps = 0; col_i = 0; 
    val_alpha = 0; alpha_idx = 0;
    Cb = setdiff((1:n),C);
    L,U,P = lu(A[R,C]);
    Arp = A[R[P],:];
    flag = true;
    while flag
        flag = false;
        if TP == :FI || TP == :FP
            for (i,Cb_i) in enumerate(Cb)
                # Get multiplicator alpha from LU Factoration
                alpha = abs.(U\(L\Arp[:,Cb_i]));
                if TP == :FI
                    alpha_idx = findfirst(x -> x >= 1+eps, alpha); 
                elseif TP == :FP
                    val_alpha,alpha_idx = findmax(alpha); 
                    if val_alpha < 1 + eps
                        alpha_idx = nothing
                    end
                end
                if !isnothing(alpha_idx)
                    col_i = i;
                    flag = true;
                    break
                end
            end
        elseif TP == :BI 
            # Get multiplicator alpha from LU Factoration
            alpha = abs.(U\(L\Arp[:,Cb]));
            val_alpha,alpha_idx = findmax(alpha); 
            if val_alpha < 1 + eps
                alpha_idx = nothing
            end
            if !isnothing(alpha_idx)
                col_i = alpha_idx[2];
                alpha_idx = alpha_idx[1];
                flag = true;
            end
        else
            error("LS type $(TP) does not exist.")
        end
        if flag
            # increase number of swaps
            swaps += 1;
            # update LU with rank-one update
            sub_cols = Arp[:,Cb[col_i]]-Arp[:,C[alpha_idx]];
            ei = zeros(r); ei[alpha_idx] = 1;
            L,U = lu_rank_one_update(L, U, sub_cols, ei, r);
            # swap elements elements from C to Cb and vice-versa
            C[alpha_idx],Cb[col_i] = Cb[col_i],C[alpha_idx]; 
        end
    end
    # get output
    time_ls = (time_ns() - time_start)/1e9;
    det_Ar = det(A[R,C]);  
    A_hat =  A[:,C];
    H_hat = (A_hat'*A_hat) \ A_hat';
    H = zeros(n,m);
    H[C,:] .= H_hat;
    admmsol = SolutionADMM();
    admmsol.time = time_ls;
    admmsol.H = H;
    admmsol.iter = swaps;
    admmsol.z = det_Ar;
    return admmsol,C
end

function LS_21(A,C,TP::Symbol)
    time_start = time_ns()
    m,n = size(A); r = length(C);
    swaps = 0;
    A_hat =  A[:,C];
    H_hat = (A_hat'*A_hat) \ A_hat';
    norm21 = getnorm21(H_hat);
    Cb = setdiff((1:n),C);
    Psi = H_hat*A[:,Cb];
    i_save,j_save = 0,0; 
    vbar_save,vj_save = Nothing,0;
    norm_Hj_save = 0.0;
    W = H_hat*H_hat'
    norm_H_hat_row = norm.(eachrow(H_hat)).^2;
    flag = true;
    while flag
        flag = false;
        for col_i in (1:n-r)
            # Get vector v to efficiently compute new 21-norm
            v = Psi[:,col_i] # H_hat*A[:,Cb_i];
            Sj = collect(2:r); # remove j from (1:r)
            for j = (1:r)
                if j >= 2
                    Sj[j-1] = j-1;
                end
                vj = v[j];
                if abs(vj) < 1e-8
                    continue;
                end
                v_bar = -v[Sj]./vj;
                norm_Hj = sqrt(norm_H_hat_row[j])
                new_21norm  = norm_Hj/abs(vj);
                new_21norm += sum(sqrt.(norm_H_hat_row[Sj] +  2*v_bar.*W[Sj,j] +  v_bar.^2*norm_Hj^2 )); 
                if new_21norm < norm21
                    flag = true; 
                    norm21 = new_21norm;
                    i_save,j_save = col_i,j;
                    vbar_save = v_bar;
                    norm_Hj_save = norm_Hj;
                    vj_save = vj;
                    if TP == :FI
                        break  
                    end              
                end
            end
            if flag && (TP == :FP || TP == :FI)
                break
            end
        end
        if flag
            # increase number of swaps
            swaps += 1;
            # update H_hat efficiently
            i,j = i_save,j_save;
            Sj = collect(1:r);
            deleteat!(Sj,j);
            H_hat[Sj,:] .=  H_hat[Sj,:] .+ vbar_save.*H_hat[j,:]';
            H_hat[j,:] .= (1/vj_save)*H_hat[j,:];
            # update columns
            C[j],Cb[i] = Cb[i],C[j]; 
            # update W
            Y = W[Sj,j]*vbar_save';
            W[Sj,Sj] += (vbar_save*vbar_save')*norm_Hj_save^2 + Y + Y';
            W[Sj,j] = (1/vj_save)*(W[Sj,j] + (norm_Hj_save^2)*vbar_save)
            W[j,Sj] =  W[Sj,j];
            W[j,j] = (norm_Hj_save/vj_save)^2
            # update Psi
            Psi[Sj,:] += vbar_save.*Psi[j,:]';
            Psi[j,:] *= (1/vj_save);
            Psi[:,i] .= H_hat*A[:,Cb[i]]; 
            # update norm H hat squared
            norm_H_hat_row = norm.(eachrow(H_hat)).^2;
        end
    end
    # get output
    time_ls21 = (time_ns() - time_start)/1e9;
    A_hat =  A[:,C];
    H_hat = (A_hat'*A_hat) \ A_hat';
    H = zeros(n,m);
    H[C,:] .= H_hat;
    admmsol = SolutionADMM();
    admmsol.time = time_ls21;
    admmsol.H = H;
    admmsol.iter = swaps;
    admmsol.z = getnorm21(admmsol.H);
    return admmsol
end


function LS_det_not_eff(A,R,C,TP::Symbol)
    time_start = time_ns()
    m,n = size(A); r = length(C);
    eps = 1e-4;
    swaps = 0; col_i = 0; 
    val_alpha = 0; alpha_idx = 0;
    Cb = setdiff((1:n),C);
    det_best = abs(det(A[R,C]));
    C_best = deepcopy(C)
    C2 = deepcopy(C);
    flag = true
    idx_ib,idx_jb = 0,0;
    # @show Cb
    while flag
        flag = false
        det_best2 = deepcopy(det_best)
        for (idx_i,i) in enumerate(Cb)
            for (idx_j,j) in enumerate(C)
                C2 = deepcopy(C);
                C2[idx_j] = i;
                if abs(det(A[R,C2])/det_best) >=  1 + eps
                    flag = true;
                    C_best = deepcopy(C2)
                    det_best = abs(det(A[R,C2]))
                    idx_ib,idx_jb = idx_i,idx_j
                    if TP ==:FI
                        break
                    end
                end
            end
            if (TP ==:FI || TP == :FP) && flag
                break
            end
        end
        if flag 
            # C = deepcopy(C_best)
            C[idx_jb],Cb[idx_ib] = Cb[idx_ib],C[idx_jb]; 
            swaps += 1
        end

    end
    time_ls = (time_ns() - time_start)/1e9;
    det_Ar = det(A[R,C]);  
    A_hat =  A[:,C];
    H_hat = (A_hat'*A_hat) \ A_hat';
    H = zeros(n,m);
    H[C,:] .= H_hat;
    admmsol = SolutionADMM();
    admmsol.time = time_ls;
    admmsol.H = H;
    admmsol.iter = swaps;
    admmsol.z = det_Ar;
    return admmsol,C
end


function LS_21_not_eff(A,C,TP::Symbol)
    time_start = time_ns()
    m,n = size(A); r = length(C);
    eps = 1e-4;
    swaps = 0; col_i = 0; 
    val_alpha = 0; alpha_idx = 0;
    Cb = setdiff((1:n),C);
    n21_best = getnorm21(pinv(A[:,C]));
    C_best = deepcopy(C)
    C2 = deepcopy(C);
    flag = true
    idx_ib,idx_jb = 0,0;
    # @show Cb
    while flag
        flag = false
        # det_best2 = deepcopy(det_best)
        for (idx_i,i) in enumerate(Cb)
            for (idx_j,j) in enumerate(C)
                C2 = deepcopy(C);
                C2[idx_j] = i;
                if getnorm21(pinv(A[:,C2])) < n21_best #- eps 
                # if abs(det(A[R,C2])/det_best) >=  1 + eps
                    flag = true;
                    C_best = deepcopy(C2)
                    n21_best = getnorm21(pinv(A[:,C2]))#abs(det(A[R,C2]))
                    idx_ib,idx_jb = idx_i,idx_j
                    if TP ==:FI
                        break
                    end
                end
            end
            if (TP ==:FI || TP == :FP) && flag
                break
            end
        end
        if flag 
            # C = deepcopy(C_best)
            C[idx_jb],Cb[idx_ib] = Cb[idx_ib],C[idx_jb]; 
            swaps += 1
        end

    end
    # get output
    time_ls21 = (time_ns() - time_start)/1e9;
    A_hat =  A[:,C];
    H_hat = (A_hat'*A_hat) \ A_hat';
    H = zeros(n,m);
    H[C,:] .= H_hat;
    admmsol = SolutionADMM();
    admmsol.time = time_ls21;
    admmsol.H = H;
    admmsol.iter = swaps;
    admmsol.z = getnorm21(admmsol.H);
    return admmsol
end