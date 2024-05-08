using Convex, LinearAlgebra, SCS,Gurobi,ForwardDiff

m,n,r = 5,4,3;

A = rand(m,r)*rand(r,n);
Y = rand(r,m) - 2*rand(r,m);
W = rand(n,m-r) - 2*rand(n,m-r);
rho = 5*rand();
rho = 1;
lambda = 1/(rho);
# lambda = q.λ[end]
U,S,V = svd(A,full=true);
U1 = U[:,1:r]; U2 = U[:,r+1:end];
V1 = V[:,1:r]; V2 = V[:,r+1:end];
E = Variable(n,m)
expr = lambda*norm(vec(E),1) + (1/2)*sumsquares(vec(V1'*E-Y)) + (1/2)*sumsquares(vec(W - E*U2));
problem = minimize(expr)
solve!(problem, Gurobi.Optimizer)#; silent_solver=true)
evaluate(expr)
fe(E) =  lambda*norm(E,1) + (1/2)*norm(V1'*E-Y)^2 + (1/2)*norm(W - E*U2)^2
E = E.value
fe(E)


# E2 = closed_form1n(0.5*V1*Y*U2*U2' + 0.5*V1*V1'*W*U2',lambda)
# norm(E-E2)

# 
# f2e(E) = (1/4)*norm(V1'*E-Y)^2 + (1/4)*norm(E*U2 - W)^2
# ForwardDiff.gradient(fe,E2)
# X =0.5*V1*(Y*U2 + V1'*W)*U2';



# @show fe(E),fe(E2)

# my_grad(E) = lambda*sign.(E) + 0.5*V1*V1'*E - 0.5*V1*Y + 0.5*E*U2*U2' - 0.5*W*U2';
# my_grad(E) = lambda*V1'*sign.(E) + 0.5*V1'*E - 0.5*Y + 0.5*V1'*E*U2*U2' - 0.5*V1'*W*U2';
# my_grad(E) = lambda*V1'*sign.(E)*U2 + 0.5*V1'*E*U2 - 0.5*Y*U2 + 0.5*V1'*E*U2 - 0.5*V1'*W;
# my_grad(E) = lambda*V1'*sign.(E)*U2 + V1'*E*U2 - 0.5*Y*U2 - 0.5*V1'*W;
# my_grad(E) = V1'*(lambda*sign.(E)+E)*U2 - 0.5*(Y*U2 + V1'*W);

# E2 = X - lambda*sign.(X)
# for i = (1:n)
#     for j = (1:m)
#         if abs(X[i,j]) <= lambda
#             E2[i,j] = 0.0;
#         else
#             E2[i,j] = X[i,j] - lambda*sign(X[i,j]);
#         end
#     end
# end

# E2 = closed_form1n(X,lambda)
# my_grad2(E) = 0.5*V1*V1'*E - 0.5*V1*Y + 0.5*E*U2*U2' - 0.5*W*U2';


# z1 = norm(vec(E),1) + (rho/2)*norm(vec(V1'*E-Y))^2 + (rho/2)*norm(vec(E*U2 - W))^2
# z2 = norm(vec(E2),1) + (rho/2)*norm(vec(V1'*E2-Y))^2 + (rho/2)*norm(vec(E2*U2 - W))^2
# z1,z2



# z1 = norm(vec(D),1) + (a/2)*norm(vec(D-A))^2 + (b/2)*norm(vec(D-B))^2

# D2 = closed_form1n(2A,1/(1))
# ab = a+b;ab_inv = 1/ab;AB = (ab_inv)*(a*A+b*B); 
# D2 = sign.(AB).*max.(abs.(AB) .- ab_inv, 0)
# z2 = norm(vec(D2),1) + (a/2)*norm(vec(D2-A))^2 + (b/2)*norm(vec(D2-B))^2
# (z1-z2)