using LinearAlgebra, Gurobi, MosekTools, JuMP

function testupdate210(n,r,V1Dinv,V2,Y)
    # n,r = inst.n,inst.r;
    model = Model(Mosek.Optimizer);
    # set_time_limit_sec(model, 15)
    set_silent(model);
    # set_attribute(model, "TimeLimit", 15)
    # set_attribute(model, "Presolve", 0)
    # set_attribute(model, "BarConvTol", 1e-5)
    # set_attribute(model, "BarQCPConvTol", 1e-5)
    # set_attribute(model, "FeasibilityTol", 1e-5)
    # set_attribute(model, "OptimalityTol", 1e-5)
    # set_silent(model);
    @variable(model, t[i in (1:n)]);
    # @variable(model, Z[1:n-r,1:r]);
    @variable(model, W[1:n,1:r]);
    Cmat = V1Dinv + W #V2*Z
    @constraint(model, [i = (1:n)], [t[i];  Cmat[i,:]] in SecondOrderCone())
    M = Cmat - Y
    @objective(model,Min,sum(t) + 0.5*sum(M[i,j]^2 for i = (1:n), j = (1:r)))
    optimize!(model);
    # Z = value.(Z)
    W = value.(W)
    return W
end

function testupdate21(n,r,V1Dinv,V2,Y)
    # n,r = inst.n,inst.r;
    model = Model(Mosek.Optimizer);
    # set_time_limit_sec(model, 15)
    set_silent(model);
    @variable(model, t[i in (1:n)]);
    # @variable(model, Z[1:n-r,1:r]);
    @variable(model, W[1:n,1:r]);
    Cmat = W #V2*Z
    @constraint(model, [i = (1:n)], [t[i];  Cmat[i,:]] in SecondOrderCone())
    M = Cmat - Y
    @objective(model,Min,sum(t) + 0.5*sum(M[i,j]^2 for i = (1:n), j = (1:r)))
    optimize!(model);
    # Z = value.(Z)
    W = value.(W)
    return W
end


# Y = rand(n,r) - 2*rand(n,r)
# rho = 1
# rho_inv = 1
# E = zeros(n,r)
# norm_rows_y = norm.(eachrow(Y))
# svec2normW = sortperm(norm_rows_y,rev=true)
# for (i,idx) in enumerate(svec2normW)
#     el = norm_rows_y[idx]
#     if rho_inv < el
#         E[i,:] = ((el - rho_inv)/el) * Y[i,:]
#     end
# end
# Ztemp = V2'*(E-V1Dinv)
# Ztrue = testupdate210(n,r,V1Dinv,V2,Y)

# objvaltemp = getnorm21(V1Dinv + V2*Ztemp) + 0.5*norm(V1Dinv + V2*Ztemp - Y)^2
# objvaltrue = getnorm21(V1Dinv + V2*Ztrue) + 0.5*norm(V1Dinv + V2*Ztrue - Y)^2

Y = rand(n,r) - 2*rand(n,r)
rho = 1
rho_inv = 1
E = zeros(n,r)
norm_rows_y = norm.(eachrow(Y))
svec2normW = sortperm(norm_rows_y,rev=true)
# for (i,idx) in enumerate(svec2normW)
#     el = norm_rows_y[idx]
#     if rho_inv < el
#         E[i,:] = ((el - rho_inv)/el) * Y[i,:]
#     end
# end
for i = (1:n)
    yi = Y[i,:]
    el = norm(yi)
    if rho_inv < el
        E[i,:] = ((el - rho_inv)/el) * Y[i,:]
    end
end

Wtemp = (E-V1Dinv)
Wtrue = testupdate210(n,r,V1Dinv,V2,Y)

objvaltemp = getnorm21(V1Dinv + Wtemp) + 0.5*norm(V1Dinv + Wtemp - Y)^2
objvaltrue = getnorm21(V1Dinv + Wtrue) + 0.5*norm(V1Dinv + Wtrue - Y)^2

# objvaltemp = getnorm21(E) + 0.5*norm(E - Y)^2
# objvaltrue = getnorm21(Etrue) + 0.5*norm(Etrue - Y)^2

objvaltemp,objvaltrue
