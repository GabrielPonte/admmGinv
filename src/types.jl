mutable struct GinvInst
    A::Matrix{Float64}
    m::Int64
    n::Int64
    r::Int64
end

mutable struct GinvInit
    U1::Matrix{Float64}
    U2::Matrix{Float64}
    V1::Matrix{Float64}
    V2::Matrix{Float64}
    Dinv::Vector{Float64}
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

mutable struct SolutionExtendedADMM
    z::Float64
    iter::Int64
    time::Float64
    H::Matrix{Float64}
    res_cpE::Float64
    res_cdE::Float64
    res_cpB::Float64
    res_cdB::Float64
    eps_pE::Float64
    eps_dE::Float64
    eps_pB::Float64
    eps_dB::Float64

    SolutionExtendedADMM() = new()
end

mutable struct GinvResultSolver
    m::Int64
    n::Int64
    r::Int64
    objval::Float64
    NZC::Int64
    NZR::Int64
    norm_0::Int64
    norm_1::Float64
    norm_21::Float64
    iter::Int64
    time::Float64
    norm_fro::Float64
    bool_p1::Bool
    bool_p2::Bool
    bool_p3::Bool
    bool_p4::Bool
    p1::Float64
    p2::Float64
    p3::Float64
    p4::Float64

    GinvResultSolver() = new()

    GinvResultSolver(
        m,n,r,objval,
        NZC,NZR,norm_0,norm_1,norm_21,iter,time,norm_fro,
        bool_p1,bool_p2,bool_p3,bool_p4,
        p1,p2,p3,p4
    ) = new(
        m,n,r,objval,
        NZC,NZR,norm_0,norm_1,norm_21,iter,time,norm_fro,
        bool_p1,bool_p2,bool_p3,bool_p4,
        p1,p2,p3,p4
    )
end

mutable struct GinvResultADMM
    m::Int64
    n::Int64
    r::Int64
    objval::Float64
    NZC::Int64
    NZR::Int64
    norm_0::Int64
    norm_1::Float64
    norm_21::Float64
    iter::Int64
    time::Float64
    norm_fro::Float64
    res_pri::Float64
    res_dual::Float64
    res_d_ML::Float64
    res_opt::Float64
    eps_p::Float64
    eps_d::Float64
    bool_p1::Bool
    bool_p2::Bool
    bool_p3::Bool
    bool_p4::Bool
    p1::Float64
    p2::Float64
    p3::Float64
    p4::Float64

    GinvResultADMM() = new()

    GinvResultADMM(
        m,n,r,objval,
        NZC,NZR,norm_0,norm_1,norm_21,iter,time,norm_fro,
        res_pri,res_dual,res_d_ML,res_opt,eps_p,eps_d,
        bool_p1,bool_p2,bool_p3,bool_p4,
        p1,p2,p3,p4
    ) = new(
        m,n,r,objval,
        NZC,NZR,norm_0,norm_1,norm_21,iter,time,norm_fro,
        res_pri,res_dual,res_d_ML,res_opt,eps_p,eps_d,
        bool_p1,bool_p2,bool_p3,bool_p4,
        p1,p2,p3,p4
    )
end

mutable struct GinvResultExtendedADMM
    m::Int64
    n::Int64
    r::Int64
    objval::Float64
    NZC::Int64
    NZR::Int64
    norm_0::Int64
    norm_1::Float64
    norm_21::Float64
    iter::Int64
    time::Float64
    norm_fro::Float64
    res_cpE::Float64
    res_cdE::Float64
    res_cpB::Float64
    res_cdB::Float64
    eps_pE::Float64
    eps_dE::Float64
    eps_pB::Float64
    eps_dB::Float64
    bool_p1::Bool
    bool_p2::Bool
    bool_p3::Bool
    bool_p4::Bool
    p1::Float64
    p2::Float64
    p3::Float64
    p4::Float64

    GinvResultExtendedADMM() = new()

    GinvResultExtendedADMM(
        m,n,r,objval,
        NZC,NZR,norm_0,norm_1,norm_21,iter,time,norm_fro,
        res_cpE,res_cdE,res_cpB,res_cdB,eps_pE,eps_dE,eps_pB,eps_dB,
        bool_p1,bool_p2,bool_p3,bool_p4,
        p1,p2,p3,p4
    ) = new(
        m,n,r,objval,
        NZC,NZR,norm_0,norm_1,norm_21,iter,time,norm_fro,
        res_cpE,res_cdE,res_cpB,res_cdB,eps_pE,eps_dE,eps_pB,eps_dB,
        bool_p1,bool_p2,bool_p3,bool_p4,
        p1,p2,p3,p4
    )
end