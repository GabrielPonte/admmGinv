using Lasso, LinearAlgebra
H = [kron(Matrix(I,m,m),V1');kron(U2',Matrix(I,n,n))];
g = [vec(Y);vec(W)];
lambda = 0.1*rand()
lista_lamb = [lambda/size(H,1)]#collect(1000:-1:1)
q = fit(LassoPath,H, g,
        intercept = false,
        standardize=false,
        # dofit = false,
        # algorithm=NaiveCoordinateDescent,
        λ=lista_lamb,
        α=1.0)
x = q.coefs[:,end]

fl(x) = lambda*norm(x,1) + (1/(2))*norm(H*x-g)^2
v = Variable(n*m)
expr = lambda*norm(v,1) + (1/(2))*sumsquares(H*v-g);
problem = minimize(expr)
solve!(problem, Gurobi.Optimizer; silent_solver=true)
v = v.value
fl(x),fl(v)
# E2 = reshape(q.coefs[:,end],n,m);
# fe(E)

