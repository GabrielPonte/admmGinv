using Convex, LinearAlgebra, SCS,Gurobi

A = 2*(rand(3,3)-2*rand(3,3));
B = 2*(rand(3,3)-2*rand(3,3));
D = Variable(3,3);
a,b = 2*rand(),2*rand()
# a,b=3,3
expr = norm(vec(D),1) + (a/2)*sumsquares(vec(D-A)) + (b/2)*sumsquares(vec(D-B));
problem = minimize(expr)
solve!(problem, Gurobi.Optimizer; silent_solver=true)
evaluate(expr)
D = D.value

z1 = norm(vec(D),1) + (a/2)*norm(vec(D-A))^2 + (b/2)*norm(vec(D-B))^2

# D2 = closed_form1n(2A,1/(1))
ab = a+b;ab_inv = 1/ab;AB = (ab_inv)*(a*A+b*B); 
D2 = sign.(AB).*max.(abs.(AB) .- ab_inv, 0)
z2 = norm(vec(D2),1) + (a/2)*norm(vec(D2-A))^2 + (b/2)*norm(vec(D2-B))^2
(z1-z2)