function [C,time] = qrGL(A,m,n,r,minSigmaFactor)
    tic
    % Start Greedy Light 
    C = [];
    sizeC = length(C);
    Qc = eye(m,m);
    Rc = zeros(m,r);
    while sizeC < r
        flag = 0;
        v = zeros(r,1);
        v(sizeC+1) = 1;
        for i = (1:n)
            X = isempty(find(C==i));
            if X==1
                [Qi,Ri] = qrupdate(Qc,Rc,A(:,i),v);
                min_sigma = abs(Ri(sizeC+1,sizeC+1));
                % C_save = [C;i];
                % min_sigma = min(svd(A(:,C_save),'econ'));
                if min_sigma > 10^(-minSigmaFactor)
                    flag = 1;
                    C = [C;i];
                    Qc = Qi;
                    Rc = Ri;
                    break
                end
            end
        end
        if flag == 0
           minSigmaFactor = minSigmaFactor + 1; 
        end
        sizeC = length(C);
    end % end Greedy Light
    time = toc;
end