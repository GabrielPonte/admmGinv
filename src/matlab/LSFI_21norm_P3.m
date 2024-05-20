function [norm_val,elapsed_time,swaps,C_out] = LSFI_21norm_P3(A,r,n,C)

    tic
    
    m = size(A,1);
    Cb = setdiff((1:n),C); % complement of C
    A_hat = A(:,C);
    H_hat = (A_hat'*A_hat) \ A_hat'; % H_hat = pinv(A_hat);
    value_ref = sum(sqrt(sum(H_hat.^2, 2)));
    
    swaps = 0;
    
    flag = true;
    
    %% LOCAL SEARCH 2,1-norm

    while flag
        
        flag = false;
        
        for column = (1:n-r)
            
            v = H_hat*A(:,Cb(column));
            
            for j = (1:r)
                H_hat = pinv(A(:,C));
                Sj = [(1:j-1) (j+1:r)]; % remove j from (1:r)
                vj = v(j);
                
                if abs(vj) < 1e-8
                    continue;
                end
                
                v_bar = -v(Sj)/vj;
                Hs = H_hat(Sj,:);
                Hj = H_hat(j,:);
                norm_Hj = norm(Hj,2);
                w = Hs*Hj';
                new_21norm = norm_Hj/abs(vj) +  sum(sqrt(sum(Hs.^2, 2) +  2*v_bar.*w +  (v_bar).^2*norm_Hj^2 ));

                if new_21norm < value_ref
                    flag = true; 
                    swaps = swaps + 1;
                    % update H_hat efficiently
                    H_hat(Sj,:) = Hs + v_bar.*Hj;
                    H_hat(j,:) = (1/vj)*Hj;
                    % update columns
                    Cb_save = Cb(column);
                    Cb(column) = C(j); 
                    C(j) = Cb_save;
                    value_ref = new_21norm;
                    break                
                end
                
            end
            if flag
                break
            end
        end
    end
    
    % A_hat = A(:,C);
    % H_hat = inv(A_hat'*A_hat) * A_hat';
    real_21norm = sum(sqrt(sum(pinv(A(:,C)).^2, 2)));
    abs(real_21norm-value_ref)
    elapsed_time = toc;
    C_out = C;
    norm_val = value_ref;
    
end