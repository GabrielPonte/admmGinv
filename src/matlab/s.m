function [norm,time,swaps,C_out] = LSFI_Det_P3(A,r,n,R,C)

    tic;
    
    swaps = 0;
    
    %% Knowing the rows and columns that i'm not using.
    %  The current rows = R and the current columns = C
    
    Cb = [];
    
    Ar = A(R,C);
    
    for j = (1:n)
        
        n1 = find(C==j);
        
        if isempty(n1)==1
            
            Cb = [Cb;j];
        end
    end
    
    flag = 1;
    
    while flag > 0
        
        flag = 0;
        
        %% Local Search
        
        % FOR COLUMNS
        
        Ar = A(R,C);
    
        [L2,U2,P2] = lu(Ar);
        
        for i = (1: n-r)
            
            % LU Factoration
            
            b = P2 * A(R,Cb(i));
            
            y = L2\b;
            
            alfa = U2\y;
            
            
            % Changing Ar
            
            local_alfa = find(abs(alfa)>1);
            
            T = isempty(local_alfa);
            
            if T==0
                
                local_alfa = local_alfa(1);
                
                swaps = swaps + 1;
                
                el_save = C(local_alfa) ;
                
                C(local_alfa) = Cb(i);
                
                Cb(i) = el_save;
                
                Ar = A(R,C);
                
                [L2,U2,P2] = lu(Ar);
                
                flag = flag+1;
            end
            
        end
    end
    
    time = toc;
    A_hat =  A(:,C);
    
    H_hat = (A_hat'*A_hat) \ A_hat';
    % H_hat = pinv(A_hat);
    
    norm = sum(sum(abs(H_hat)));
    C_out = C;
    
end