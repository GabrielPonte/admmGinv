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

function subprobD(Y,rho)
    n,r = size(Y);
    model = Model(Mosek.Optimizer)
    set_silent(model)
    @variable(model,E[1:n,1:r])
    @variable(model,Z[1:n,1:m])
    @variable(model,D[1:n,1:m])
    EU = [E zeros(n,m-r)]*U';
    for i = (1:n)
        for j =(1:m)
            @constraint(model,Z[i,j] >= D[i,j])
            @constraint(model,Z[i,j] >= -D[i,j])
            @constraint(model,D[i,j] == EU[i,j])
            
        end
    end
    
    
    M = D-[Y zeros(n,m-r)]*U'
    M2 = D-Y*U1';
    @objective(model,Min, sum(Z[i,j] for i = (1:n), j =(1:m)) + (rho/2)*sum(M[i,j]^2 for i = (1:n), j =(1:m)))
    optimize!(model)
    @show objective_value(model)
    @show norm(value.(D))
    D=value.(D)
    E = value.(E)
    # E = value.(E)
    # X = value.(X)
    # @show norm(E - X*U1)
    # @show rank(X)
    # W = V2*Z*U1'
    return 

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


rho = 1e-4
m,n,r = 21,15,10;
# A = rand((1:3),m,r)*rand((1:4),r,n);
A = rand(m,r)*rand(r,n) - 2*rand(m,r)*rand(r,n);
U = svd(A,full=true).U;
U1,U2 = U[:,1:r],U[:,r+1:end];
U,Σ,V = svd(A,full=true);
D = Σ[1:r];
Dinv = 1 ./ D;
U1 = U[:,1:r];
V1,V2 = V[:,1:r],V[:,r+1:n];
Y = 199*(rand(n,r) - 2*rand(n,r));
rho = 1
subprobE1(Y,rho);
#subprobD(Y,rho);

# 1

