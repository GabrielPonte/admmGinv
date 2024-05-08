using Gurobi, JuMP, LinearAlgebra

# D = rand(5,5);
# z = rand(5);
model = Model(Gurobi.Optimizer)
@variable(model,x[1:5]);
@variable(model,y[1:5]);
@constraint(model,y .>= D*x);
@constraint(model,y .>=-D*x);
@objective(model,Min,sum(y) + 0.5*sum((x[i]-z[i])^2 for i =(1:5)))
optimize!(model)
z1 = objective_value(model)


model = Model(Gurobi.Optimizer)
@variable(model,u[1:5]);
# @variable(model,y[1:5]);
# @constraint(model,y .>= D*x);
# @constraint(model,y .>=-D*x);
Du = D'*u;
@objective(model,Min,0.5*sum((-Du[i] + z[i])^2 for i =(1:5)))
@constraint(model,u .<=  1);
@constraint(model,u .>= - 1);
optimize!(model)
z2 = objective_value(model);

z1,z2