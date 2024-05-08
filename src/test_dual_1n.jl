using JuMP, Gurobi, LinearAlgebra, MosekTools

function toU(Y,rho)
    n,r = size(Y);
    model = Model(Mosek.Optimizer)
    set_silent(model)
    @variable(model,E[1:n,1:r])
    @variable(model,X[1:n,1:r])
    EU = E#*U1';
    for i = (1:n)
        for j =(1:r)
            @constraint(model,X[i,j] >= EU[i,j])
            @constraint(model,X[i,j] >= -EU[i,j])
        end
    end
    M = E-Y
    @objective(model,Min, sum(X[i,j] for i = (1:n), j =(1:r)) + (rho/2)*sum(M[i,j]^2 for i = (1:n), j =(1:r)))
    optimize!(model)
    @show objective_value(model)
    E = value.(E)
    # W = V2*Z*U1'
    return E
end

function subprobE1(Y,rho)
    n,r = size(Y);
    model = Model(Mosek.Optimizer)
    set_silent(model)
    @variable(model,E[1:n,1:r])
    @variable(model,Z[1:n,1:m])
    @variable(model,X[1:n,1:m])
    EU = E*U1';
    for i = (1:n)
        for j =(1:m)
            @constraint(model,Z[i,j] >= X[i,j])
            @constraint(model,Z[i,j] >= -X[i,j])
           @constraint(model,X[i,j] == EU[i,j])
        end
    end
    
    M = X-Y*U1'
    @objective(model,Min, sum(Z[i,j] for i = (1:n), j =(1:m)) + (rho/2)*sum(M[i,j]^2 for i = (1:n), j =(1:m)))
    optimize!(model)
    @show objective_value(model)
    E = value.(E)
    X = value.(X)
    @show norm(E - X*U1)
    @show rank(X)
    # W = V2*Z*U1'
    return X
end

function subprobD(Y,U1,rho)
    m = size(U1,1);
    n,r = size(Y);
    model = Model(Mosek.Optimizer)
    set_silent(model)
    @variable(model,E[1:n,1:r])
    @variable(model,Z[1:n,1:m])
    # @variable(model,D[1:n,1:m])
    EU = E*U1';
    # DU2 = D*U2;
    # @constraint(model,D*U2 .== 0)
    for i = (1:n)
        for j =(1:m)
            @constraint(model,Z[i,j] >= EU[i,j]);
            @constraint(model,Z[i,j] >= -EU[i,j]);            
        end
    end
    
    
    M = E - Y;#D-Y*U1'
    # M2 = D-Y*U1';
    @objective(model,Min, sum(Z[i,j] for i = (1:n), j =(1:m)) + (rho/2)*sum(M[i,j]^2 for i = (1:n), j =(1:r)))
    optimize!(model)
    
    # @show norm(value.(D))
    # D = value.(D)
    # E = D*U1;
    E = value.(E)
    #@show abs(objective_value(model) - (norm(E*U1',1) + (rho/2)*norm(E - Y)^2))
    #D = E*U1';
    # X = value.(X)
    # @show norm(E - X*U1)
    # @show rank(X)
    # W = V2*Z*U1'
    return E#,D

end

function subprobE2(Y,rho)
    n,r = size(Y);
    model = Model(Mosek.Optimizer)
    set_silent(model)
    @variable(model,E[1:n,1:r])
    @variable(model,Z[1:n,1:m])
    @variable(model,X[1:n,1:m])
    EU = [E zeros(n,m-r)]*U';
    for i = (1:n)
        for j =(1:m)
            @constraint(model,Z[i,j] >= X[i,j])
            @constraint(model,Z[i,j] >= -X[i,j])
           @constraint(model,X[i,j] == EU[i,j])
        end
    end
    
    M = X-Y*U1'
    @objective(model,Min, sum(Z[i,j] for i = (1:n), j =(1:m)) + (rho/2)*sum(M[i,j]^2 for i = (1:n), j =(1:m)))
    optimize!(model)
    @show objective_value(model)
    E = value.(E);
    X = value.(X);
    @show norm(E)
    @show rank(X)
    # W = V2*Z*U1'
    return 
end

using Convex, LinearAlgebra, SCS,Gurobi

function ALM_subprobE(rho,eta,Y,U1,U2)
    Phi = ones(n,m-r);
    D = Nothing
    G = Y*U1';
    k = 0;
    while true
        re = rho+eta;re_inv = 1/re;AB = (re_inv)*(rho*G+eta*Phi*U2'); 
        # if iter <= 35
        D = (sign.(AB).*max.(abs.(AB) .- re_inv, 0));
        # else
        #     D = Variable(n,m);
        #     expr = norm(vec(D),1) + (rho/2)*sumsquares(vec(D-G)) + (eta/2)*sumsquares(vec(D-Phi*U2'));
        #     problem = minimize(expr)
        #     solve!(problem, Gurobi.Optimizer; silent_solver=true)
        #     evaluate(expr)
        #     D = D.value
        # end
        DU2 = D*U2;
        Phi = Phi - DU2;
        k = k +1;
        prim_res = norm(DU2);
        # @show k,prim_res
        # @show k,norm(DU2)
        if prim_res < 1e-4
            E = D*U1;
            # D1 = copy(D);
            # z1 = norm(D,1) + (rho/2)*norm(D-G)^2 + (eta/2)*norm(D*U2-Phi)^2
            # E2313,D = subprobD(Y,U1,U2,rho);
            # z2 = norm(D,1) + (rho/2)*norm(D-G)^2 + (eta/2)*norm(D*U2-Phi)^2
            # D2 = copy(D);
            # @show iter,(z1-z2);
            # sleep(2)
            # # @show 
            return E#,Phi
        end
    end
end


# rho = 1e-4
# m,n,r = 6,5,3;
# # A = rand((1:3),m,r)*rand((1:4),r,n);
# A = rand(m,r)*rand(r,n) #- 2*rand(m,r)*rand(r,n);
# U = svd(A,full=true).U;
# U1,U2 = U[:,1:r],U[:,r+1:end];
# rho = 1; eta = 1;

# # U,Σ,V = svd(A,full=true);
# # D = Σ[1:r];
# # Dinv = 1 ./ D;
# # U1 = U[:,1:r];
# # V1,V2 = V[:,1:r],V[:,r+1:n];
# Y = (rand(n,r) - 2*rand(n,r));
# ALM_subprobD(m,n,r,rho,eta,Y,U1,U2)
# rho = 1
# subprobE1(Y,rho);
# #subprobD(Y,rho);

# # 1

