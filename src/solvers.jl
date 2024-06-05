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

function solve1_not_eff_grb(inst::GinvInst)
    A = inst.A;
    m,n,r = inst.m,inst.n,inst.r;
    model = Model(Gurobi.Optimizer);
    #set_attribute(model, "BarConvTol", 1e-5)
    #set_attribute(model, "BarQCPConvTol", 1e-5)
    #set_attribute(model, "FeasibilityTol", 1e-5)
    #set_attribute(model, "OptimalityTol", 1e-5)
    set_silent(model);
    # set_optimizer_attribute(model, "Method", 0)
    set_time_limit_sec(model, 7200)
    @variable(model, Z[1:n,1:m]);
    @variable(model, H[1:n,1:m]);
    @constraint(model, Z .>= H);
    @constraint(model, Z .>= -H);
    AH = A*H;
    AHA = A*H*A;
    HAAp = H*A*pinv(A);
    @constraint(model, AHA - A .<= 5e-6);
    @constraint(model, HAAp - H .<= 5e-6);
    @constraint(model, AH - (AH)' .<= 5e-6);
    @constraint(model, -(AHA - A) .<= 5e-6);
    @constraint(model, -(HAAp - H) .<= 5e-6);
    @constraint(model, -(AH - (AH)') .<= 5e-6);
    @objective(model,Min,sum(Z))
    time_to_solve = @elapsed optimize!(model);
    bsol = SolutionOptimizer();
    primal_status = termination_status(model);
    @show primal_status
    if primal_status != MOI.TIME_LIMIT
        @show primal_status
        bsol.H = value.(H);
        bsol.z,bsol.time =  objective_value(model),time_to_solve;
    else
        bsol.H = zeros(n,m);
        bsol.z,bsol.time = Inf, time_to_solve;
    end
    return bsol
end

function solve21_not_eff_msk(inst::GinvInst)
    A= inst.A;
    m,n,r = inst.m,inst.n,inst.r;
    model = Model(Mosek.Optimizer);
    set_time_limit_sec(model, 7200)
    set_silent(model);
    # set_attribute(model, "BASIS_REL_TOL_S", 1e-5)
    # set_attribute(model, "INTPNT_CO_TOL_DFEAS", 1e-5)
    # set_attribute(model, "INTPNT_CO_TOL_INFEAS", 1e-5)
    # set_attribute(model, "INTPNT_CO_TOL_MU_RED", 1e-5)
    # set_attribute(model, "INTPNT_CO_TOL_PFEAS", 1e-5)
    # set_attribute(model, "INTPNT_CO_TOL_REL_GAP", 1e-5)
    # set_attribute(model, "INTPNT_TOL_DFEAS", 1e-5)
    set_silent(model);
    @variable(model, t[i in (1:n)]);
    @variable(model, H[1:n,1:m]);
    @constraint(model, [i = (1:n)], [t[i];  H[i,:]] in SecondOrderCone())
    AHA = A*H*A
    @constraint(model, AHA - A .<= 5e-6);
    @constraint(model, -(AHA - A) .<= 5e-6);
    @objective(model,Min,sum(t))
    time_to_solve = @elapsed optimize!(model);

    csol = SolutionOptimizer();
    primal_status = termination_status(model);
    println(string("Solution status Mosek: ", primal_status))
    if primal_status != MOI.TIME_LIMIT
        csol.H = value.(H);
        csol.z,csol.time =  objective_value(model),time_to_solve;
    else
        csol.H = zeros(n,m);
        csol.z,csol.time = Inf, time_to_solve;
    end
    return csol
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
    Cmat = ginvInit.Vr*diagm(ginvInit.Dinv) + ginvInit.Vs*Z;
    @constraint(model, [i = (1:n)], [t[i];  Cmat[i,:]] in SecondOrderCone())
    @objective(model,Min,sum(t))
    time_to_solve = @elapsed optimize!(model);

    csol = SolutionOptimizer();
    primal_status = termination_status(model);
    println(string("Solution status Mosek: ", primal_status))
    if primal_status != MOI.TIME_LIMIT
        csol.H = ginvInit.G + ginvInit.Vs*value.(Z)*ginvInit.Ur';
        csol.z,csol.time =  objective_value(model),time_to_solve;
    else
        csol.H = zeros(n,m);
        csol.z,csol.time = Inf, time_to_solve;
    end
    return csol
end
