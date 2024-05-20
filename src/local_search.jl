using LinearAlgebra

include("types.jl")
include("util.jl")

function lu_rank_one_update(L, U, y, z)
    # only works for P = (1:n) otherwise y = y[P]
    # z = ei
    # y = A(R,C(i)) - A(R,C(alpha_idx))
    #y = y[P]; 
    n = length(y);
    for i = (1:n)
        U[i, i] = U[i, i] + y[i]*z[i];
        z[i] = z[i] / U[i, i]; 
        for j = i+1:n
            y[j] = y[j] - y[i]*L[j, i];
            L[j, i] = L[j, i] + z[i]*y[j];
            U[i, j] = U[i, j] + y[i]*z[j];
            z[j] = z[j] - z[i]*U[i, j];
        end
    end
    return L,U
end


function LS_det(A,R,C,TP::Symbol)
    time_start = time_ns()
    m,n = size(A); r = length(C);
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
                    alpha_idx = findfirst(x -> x >= 1+1e-8, alpha); 
                elseif TP == :FP
                    val_alpha,alpha_idx = findmax(alpha); 
                    if val_alpha < 1 + 1e-8
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
            if val_alpha < 1 + 1e-8
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
            L,U = lu_rank_one_update(L, U, sub_cols, ei);
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
    return H,time_ls,swaps,det_Ar
end

function LS_21(A,C,TP::Symbol)
    time_start = time_ns()
    m,n = size(A); r = length(C);
    swaps = 0;
    H_hat = pinv(A[:,C]);
    norm21 = getnorm21(H_hat);
    Cb = setdiff((1:n),C);
    i_save,j_save = 0,0; 
    vbar_save,vj_save = Nothing,0;
    flag = true;
    while flag
        flag = false;
        for (col_i,Cb_i) in enumerate(Cb)
            # Get vector v to efficiently compute new 21-norm
            v = H_hat*A[:,Cb_i];
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
                Hs = H_hat[Sj,:];
                Hj = H_hat[j,:];
                norm_Hj = norm(Hj,2);
                w = Hs*Hj;
                new_21norm  = norm_Hj/abs(vj);
                new_21norm += sum(sqrt.(norm.(eachrow(Hs)).^2 +  2*v_bar.*w +  v_bar.^2*norm_Hj^2 )); 
                if new_21norm < norm21
                    flag = true; 
                    norm21 = new_21norm;
                    i_save,j_save = col_i,j;
                    vbar_save = v_bar;
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
            # v1 = isapprox(abs(norm21-getnorm21(pinv(A[:,C]))),0;atol=1e-6)
            # v2= isapprox(norm(pinv(A[:,C])-H_hat),0;atol=1e-6)
            # @show v1,v2
        end
    end
    # get output
    time_ls = (time_ns() - time_start)/1e9;
    det_Ar = det(A[R,C]);  
    A_hat =  A[:,C];
    H_hat = (A_hat'*A_hat) \ A_hat';
    H = zeros(n,m);
    H[C,:] .= H_hat;
    return H,time_ls,swaps,det_Ar
end

m,n,r = 30,20,5;
A = rand(m,r)*rand(r,n);
R = (1:r); C = collect(1:r);
H,time_ls,swaps,det_FI = LS_det(A,R,C,:FI);
H,time_ls,swaps,det_FP = LS_det(A,R,C,:FP);
H,time_ls,swaps,det_BI = LS_det(A,R,C,:BI);
LS_21(A,C,:BI);
1;
# det_FI,det_FP,det_BI