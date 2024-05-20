using LinearAlgebra

include("types.jl")
include("util.jl")

function lu_rank_one_update(L, U, y, z)
    # only works for P = (1:n)
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


function LS_Det(A,R,C,TP::Symbol)
    time_start = time_ns()
    m,n = size(A); r = length(C);
    swaps = 0; col_i = 0; alpha_idx = 0;
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
            alpha = abs.(U\(L\Arp));
            val_alpha,alpha_idx = findmax(alpha); 
            if val_alpha < 1 + 1e-8
                alpha_idx = nothing
            end
            if !isnothing(alpha_idx)
                col_i = alpha_idx[2];
                alpha_idx = alpha_idx[1];
                flag = true;
                break
            end
        end

        if flag
            # increase number of swaps
            swaps += 1;
            # update LU with rank-one update
            sub_cols = Arp[:,Cb_i]-Arp[:,C[alpha_idx]];
            ei = zeros(r); ei[alpha_idx] = 1;
            L,U = lu_rank_one_update(L, U, sub_cols, ei);
            # swap elements elements from C to Cb and vice-versa
            # val_old_det = det(A[R,C])
            C[alpha_idx],Cb[i] = Cb[i],C[alpha_idx];
            # val_new_det = det(A[R,C])
            # @show val_alpha, val_new_det/val_old_det
            # @show norm(Arp[:,C] - L*U)        
        end
    end
    
    elapsed_time = (time_ns() - time_start)/1e9;
    A_hat =  A[:,C];
    H_hat = (A_hat'*A_hat) \ A_hat';
    H = zeros(n,m);
    H[C,:] .= H_hat;

    ls_21n = getnorm21(H);

    return C,ls_21n
end

m,n,r = 30,25,10;
A = rand(m,r)*rand(r,n);
R = (1:r); C = collect(1:r);
FP_Det(A,r,n,R,C)