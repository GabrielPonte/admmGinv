function solve1grb(inst::GinvInst,ginvInit::GinvInit)
    m,n,r = inst.m,inst.n,inst.r;
    model = Model(Gurobi.Optimizer);
    set_silent(model);
    set_time_limit_sec(model, 7200)
    set_attribute(model, "BarConvTol", 1e-5)
    set_attribute(model, "BarQCPConvTol", 1e-5)
    set_attribute(model, "FeasibilityTol", 1e-5)
    set_attribute(model, "OptimalityTol", 1e-5)
    @variable(model, t);
    @variable(model, Z[1:n-r,1:r]);
    F = ginvInit.V1DinvU1T + ginvInit.V2*Z*ginvInit.U1';
    @constraint(model, [t; vec(F)] in MOI.NormOneCone(1 + length(vec(F))));
    @objective(model,Min,t);
    time_to_solve = @elapsed optimize!(model);
    bsol = SolutionOptimizer();
    primal_status = termination_status(model);
    if primal_status != MOI.TIME_LIMIT
        bsol.H = ginvInit.V1DinvU1T + ginvInit.V2*value.(Z)*ginvInit.U1';
        bsol.z,bsol.time =  objective_value(model),time_to_solve;
    else
        bsol.H = zeros(n,m);
        bsol.z,bsol.time = Inf, time_to_solve;
    end
    return bsol
end

function solve21msk(inst::GinvInst,ginvInit::GinvInit)
    n,r = inst.n,inst.r;
    model = Model(Mosek.Optimizer);
    set_time_limit_sec(model, 7200)
    set_silent(model);
    set_attribute(model, "BASIS_REL_TOL_S", 1e-5)
    set_attribute(model, "INTPNT_CO_TOL_DFEAS", 1e-5)
    set_attribute(model, "INTPNT_CO_TOL_INFEAS", 1e-5)
    set_attribute(model, "INTPNT_CO_TOL_MU_RED", 1e-5)
    set_attribute(model, "INTPNT_CO_TOL_PFEAS", 1e-5)
    set_attribute(model, "INTPNT_CO_TOL_REL_GAP", 1e-5)
    set_attribute(model, "INTPNT_TOL_DFEAS", 1e-5)
    set_silent(model);
    @variable(model, t[i in (1:n)]);
    @variable(model, Z[1:n-r,1:r]);
    Cmat = ginvInit.V1Dinv + ginvInit.V2*Z;
    @constraint(model, [i = (1:n)], [t[i];  Cmat[i,:]] in SecondOrderCone())
    @objective(model,Min,sum(t))
    time_to_solve = @elapsed optimize!(model);
    csol = SolutionOptimizer();
    primal_status = termination_status(model);
    if primal_status != MOI.TIME_LIMIT
        csol.H = ginvInit.V1DinvU1T + ginvInit.V2*value.(Z)*ginvInit.U1';
        csol.z,csol.time =  objective_value(model),time_to_solve;
    else
        csol.H = zeros(n,m);
        csol.z,csol.time = Inf, time_to_solve;
    end
    return csol
end
