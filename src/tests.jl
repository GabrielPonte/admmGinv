function min1normFro(Y)
    model = Model(Gurobi.Optimizer);
    set_silent(model);
    @variable(model, E[1:n,1:m]);
    @variable(model, F[1:n,1:m]);
    @constraint(model, F .>= E);
    @constraint(model, F .>= -E);
    YE = Y-E;
    @objective(model,Min,sum(F) + (mu/2)*sum(YE[i,j]^2 for i =(1:n), j = (1:m)))
    optimize!(model);
    E = value.(E)
    return E
end

function minFroNorm(J,Q,V2)
    model = Model(Gurobi.Optimizer);
    set_silent(model);
    @variable(model, Z[1:n-r,1:r]);
    JZ = J - V2*Z;
    QZ = Q - V2*Z;
    @objective(model,Min,sum(JZ[i,j]^2 for i =(1:n), j = (1:r)) + sum(QZ[i,j]^2 for i =(1:n), j = (1:r)))
    optimize!(model);
    Z = value.(Z)
    return Z
end