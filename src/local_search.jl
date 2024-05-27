getLICols(A,r) = (qr(A, ColumnNorm()).p)[1:r]
getLIRows(A,r) = (qr(A', ColumnNorm()).p)[1:r]

function get_LI_RC(A,r)
    time_start = time_ns() 
    C = getLICols(A,r);
    time_start2 = time_ns()
    time_2 =  (time_start2 - time_start)/1e9
    println("Found columns in $(round(time_2,digits=3)) sec.")
    R = getLIRows(A[:,C],r);
    time_3 =  (time_ns() - time_start2)/1e9
    println("Found rows in $(round(time_3,digits=3)) sec.")
    elapsed_time = (time_ns() - time_start)/1e9
    return R,C,elapsed_time
end
    
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
    return H,time_ls,swaps,det_Ar,C
end

function LS_21(A,C,TP::Symbol)
    time_start = time_ns()
    m,n = size(A); r = length(C);
    swaps = 0;
    # H_hat = pinv(A[:,C]);
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
            v = Psi[:,col_i] #H_hat*A[:,Cb_i];
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
                # Hs = H_hat[Sj,:];
                # Hj = (H_hat[j,:]);
                
                # w = Hs*Hj; # W[Sj,j]
                # w =  W[Sj,j];
                # @show norm(w - W[Sj,j])

                norm_Hj = sqrt(norm_H_hat_row[j]) #norm(H_hat[j,:],2);
                
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
            # @show norm(W - H_hat*H_hat')
            # update Psi
            # @show size(Psi[Sj,:])

            Psi[Sj,:] += vbar_save.*Psi[j,:]';
            Psi[j,:] *= (1/vj_save);
            # @show size(A[:,Cb[i]])
            # @show size(H_hat*A[:,Cb[i]])
            Psi[:,i] .= H_hat*A[:,Cb[i]]; 
            # update norm H hat squared
            norm_H_hat_row = norm.(eachrow(H_hat)).^2;
            # @show norm(Psi -  H_hat*A[:,Cb])
            # display("text/plain", W- H_hat*H_hat')
            # return
        end
    end
    # get output
    time_ls21 = (time_ns() - time_start)/1e9;
    # det_Ar = det(A[R,C]);  
    A_hat =  A[:,C];
    H_hat = (A_hat'*A_hat) \ A_hat';
    H = zeros(n,m);
    H[C,:] .= H_hat;
    # @show getnorm21(H),norm21
    # @show time_ls21
    return H,time_ls21,swaps
end

function genPlotDetvs21()
    m =2000;
    n,r = floor(Int64,1.0*m),floor(Int64,0.1*m);
    nameInst = string("A_",m,"_",n,"_",r);
    println(string("\nStarting instance: m,n,r: ",m,",",n,",",r))
    A = getMatlabInstance(nameInst,"A");
    R_init,C_init,time_RC = get_LI_RC(A,r)
    min_21n = Inf
    res = []
    for type_det in [:FI,:FP,:BI]
        H,time_det,swaps_det,det_Ar,C = LS_det(A,copy(R_init),copy(C_init),type_det)
        n21old = getnorm21(H);
        H,time_ls21,swaps = LS_21(A,C,:FI);
        tdet = time_RC+time_det;
        n21new = round(getnorm21(H),digits=3);
        if n21new < min_21n
            min_21n = n21new
        end
        t21 = tdet + time_ls21
        push!(res,[type_det,n21old,n21new,tdet,t21,0,0]);
    end
    for i = (1:3)
        res[i][6] = (res[i][2] - min_21n)/min_21n
        res[i][7] = (res[i][3] - min_21n)/min_21n
    end
    return res
end
# m,n,r = 30,20,5;
# A = rand(m,r)*rand(r,n);
# R = (1:r); C = collect(1:r);
# H,time_ls,swaps,det_FI = LS_det(A,R,C,:FI);
# H,time_ls,swaps,det_FP = LS_det(A,R,C,:FP);
# H,time_ls,swaps,det_BI = LS_det(A,R,C,:BI);
# LS_21(A,C,:BI);
# 1;
# det_FI,det_FP,det_BI