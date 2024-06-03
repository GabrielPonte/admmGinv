function [C,time] = greedy_light(A,m,n,r,C,minSigmaFactor)
    tic
    % Start Greedy Light 
    sizeC = length(C);
    while sizeC < r
        flag = 0;
        for i = (1:n)
            X = isempty(find(C==i));
            if X==1
                C_save = [C;i];
                min_sigma = min(svd(A(:,C_save),'econ'));
                if min_sigma > 10^(-minSigmaFactor)
                    flag = 1;
                    C = [C;i];
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