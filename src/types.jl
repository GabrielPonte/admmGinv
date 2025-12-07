mutable struct GinvInst
    A::Matrix{Float64}
    m::Int64
    n::Int64
    r::Int64
end

mutable struct GinvInit
    m::Int64
    n::Int64
    r::Int64
    U1::Matrix{Float64}
    U2::Matrix{Float64}
    V1::Matrix{Float64}
    V2::Matrix{Float64}
    Dinv::Vector{Float64}
    V1Dinv::Matrix{Float64}
    V1DinvU1T::Matrix{Float64}
    V2V2T::Matrix{Float64}
end

mutable struct GinvADMM1
    U1::Matrix{Float64}
    V2::Matrix{Float64}
    Θ::Matrix{Float64}
    V1DinvU1T::Matrix{Float64}
    U1U1T::Matrix{Float64}
    V2V2T::Matrix{Float64}
    TP::Symbol
end

mutable struct GinvADMM21
    U1::Matrix{Float64}
    V2::Matrix{Float64}
    Θ::Matrix{Float64}
    V1Dinv::Matrix{Float64}
    V2V2T::Matrix{Float64}
    TP::Symbol
end

mutable struct GinvADMM20
    U1::Matrix{Float64}
    V2::Matrix{Float64}
    V1Dinv::Matrix{Float64}
    V2V2T::Matrix{Float64}
    nzr21::Int64
    ω21::Float64
end

mutable struct SolutionOptimizer
    z::Float64
    time::Float64
    H::Matrix{Float64}
    Z::Matrix{Float64}

    SolutionOptimizer() = new()
end

mutable struct SolutionADMM
    z::Float64
    iter::Int64
    time::Float64
    H::Matrix{Float64}
    res_pri::Float64
    res_dual::Float64
    res_d_ML::Float64
    res_opt::Float64
    eps_p::Float64
    eps_d::Float64

    SolutionADMM() = new()
end

mutable struct GinvResultSolver
    m::Int64
    n::Int64
    r::Int64
    objval::Float64
    norm_1::Float64
    norm_0::Int64
    norm_21::Float64
    norm_20::Int64
    norm_fro::Float64
    time::Float64
    iter::Int64
    

    GinvResultSolver() = new()

    GinvResultSolver(
        m,n,r,objval,
        norm_1,norm_0,norm_21,norm_20,norm_fro,time,iter,
    ) = new(
        m,n,r,objval,
        norm_1,norm_0,norm_21,norm_20,norm_fro,time,iter,
    )
end

mutable struct GinvResultADMM
    m::Int64
    n::Int64
    r::Int64
    objval::Float64
    norm_1::Float64
    norm_0::Int64
    norm_21::Float64
    norm_20::Int64
    norm_fro::Float64
    time::Float64
    iter::Int64
    
    GinvResultADMM() = new()

    GinvResultADMM(
        m,n,r,objval,
        norm_1,norm_0,norm_21,norm_20,norm_fro,time,iter,
    ) = new(
        m,n,r,objval,
        norm_1,norm_0,norm_21,norm_20,norm_fro,time,iter,
    )
end



mutable struct OptionsADMM
    eps::Float64
    max_iter::Float64
    time_limit::Float64
end