module MinEigLoss

using Zygote: @adjoint
using LinearAlgebra
using Zygote
using Random
using ITensors

export compute_ground_state, min_eig_loss, min_eigenvalue, H̄, H̄_abs, unitary, riemannian_gradient, riemannian_update, Opt, foo, riemannian_gradient_similarity, min_eig_loss_similarity, riemannian_gradient_similarity, H̄s, H̄s_abs, special, FF, mle_sys_special, mle_sys_unitary, mle_unitary, l1_sys_special, FF, RandomMatrix


include("ff.jl")
include("random.jl")



struct SpecialLinear{T,S<:AbstractMatrix{<:T}} <: AbstractMatrix{T}
    data::S
    function SpecialLinear{T,S}(data) where {T,S<:AbstractMatrix{<:T}}
        require_one_based_indexing(data)
        new{T,S}(data)
    end
end


import Base: size

# Implement the size method for the custom type
Base.size(sl::SpecialLinear) = size(sl.data)
Base.size(sl::SpecialLinear, d :: Int) = size(sl.data, d)
# Implement other necessary methods for AbstractMatrix
Base.getindex(sl::SpecialLinear, i::Int, j::Int) = getindex(sl.data, i, j)
Base.setindex!(sl::SpecialLinear, v, i::Int, j::Int) = setindex!(sl.data, v, i, j)
Base.IndexStyle(::Type{SpecialLinear}) = IndexCartesian()




# Function to compute the ground energy and ground state
function compute_ground_state(H)
    E, V = eigen(H)
    ground_energy = minimum(E)
    ground_state = V[:, argmin(E)]
    return ground_energy, ground_state
end

function H̄(H, u)
    U = kron(u, u)
    U * H * U'
end

function H̄s(H, w)
    W = kron(w, w)
    W * H * inv(W)
end

function stochastic(A)
    B = copy(A)
    for i in 1:size(A, 1)
        for j in 1:size(A, 2)
            if i != j
                B[i, j] = -abs(A[i, j])
            end
        end
    end
    return B
end

function H̄_abs(H, u)
    res = stochastic(H̄(H, u))
    return Hermitian(res)
end

function H̄s_abs(H, u)
    res = stochastic(H̄s(H, u))
    return res
end




minimum_eigenvalue(H :: Hermitian) = minimum(eigvals(H))
minimum_eigenvalue(H :: Matrix) = H |> eigvals |> real |> minimum

function min_eig_loss(H0, u)
    u = unitary(u)
    H = H̄_abs(H0, u)
    return -minimum_eigenvalue(H)
end

function min_eig_loss_similarity(H0, w)
    w = special(w)
    H = H̄s_abs(H0, w)
    return -minimum_eigenvalue(H)
end

# Unitary function and its adjoint
function unitary(u)
    return u
end

function special(w)
    return w
end

@adjoint unitary(u) = unitary(u), c -> (riemannian_gradient(u, c),)
@adjoint special(w) = special(w), c -> (riemannian_gradient_similarity(w, c),)

function riemannian_gradient(u, euc_grad)
    rg = euc_grad * u'
    return (rg - rg') / 2
end

function riemannian_gradient_similarity(v, euc_grad)
    v_inv = inv(v)
    g = euc_grad * v_inv
    v_inv2 = v_inv'v_inv
    v_inv2 ./= tr(v_inv2)
    - v_inv2 .* tr(g) + g
end

function riemannian_update(w, rg, step)
    exp(-step * rg) * w
end


function mle_special(Ha)
    e0 = eigmin(Ha) |> abs  
    function loss_func(v)
        return MinEigLoss.min_eig_loss_similarity(Ha, v) - e0
    end
    return loss_func
end


function mle_unitary(Ha)
    e0 = eigmin(Ha) |> abs  
    function loss_func(u)
        return MinEigLoss.min_eig_loss(Ha, u) - e0
    end
    return loss_func
end


function l1_sys_special(H:: AbstractMatrix, s_loc :: Int)
    s_sys = size(H, 1)
    n_rep :: Int = log(s_loc, s_sys) |> round |> Int
    function kron_w(w)
        W = kron(w, w)
        for i in 1:n_rep-2
            W = kron(W, w)
        end
        return W
    end

    function wrapper_special(w)
        w = MinEigLoss.special(w)
        W = kron_w(w)
        H̄ = (W * H * W') .|> abs
        return sum(H̄)
    end
end

function mle_sys_special(H:: AbstractMatrix, s_loc :: Int)
    s_sys = size(H, 1)
    n_rep :: Int = log(s_loc, s_sys) |> round |> Int
    e0 = eigmin(H) |> abs
    function kron_w(w)
        W = kron(w, w)
        for i in 1:n_rep-2
            W = kron(W, w)
        end
        return W
    end

    function wrapper_special(w)
        w = MinEigLoss.special(w)
        W = kron_w(w)
        W_inv = kron_w(inv(w))
        H̄ = (W * H * W_inv) .|> abs
        return (eigvals(H̄) |> real |> maximum) - e0
    end
end


function mle_sys_unitary(H:: AbstractMatrix, s_loc :: Int)
    s_sys = size(H, 1)
    n_rep :: Int = log(s_loc, s_sys) |> round |> Int
    e0 = eigmin(H) |> abs
    function kron_w(w)
        W = kron(w, w)
        for i in 1:n_rep-2
            W = kron(W, w)
        end
        return W
    end

    function wrapper_special(w)
        w = MinEigLoss.unitary(w)
        W = kron_w(w)
        H̄ = (W * H * W') .|> abs
        return (eigvals(H̄) |> real |> maximum) - e0
    end

    return wrapper_special
end

function mle_sys_unitary_free(H:: AbstractMatrix, s_loc :: Int)
    s_sys = size(H, 1)
    n_rep :: Int = log(s_loc, s_sys) |> round |> Int
    e0 = eigmin(H) |> abs
    function kron_w(w)
        W = kron(w, w)
        for i in 1:n_rep-2
            W = kron(W, w)
        end
        return W
    end

    function wrapper_special(w)
        w = MinEigLoss.unitary(w)
        W = kron_w(w)
        H̄ = (W * H * W') .|> abs
        return (eigvals(H̄) |> real |> maximum) - e0
    end

    return wrapper_special
end

module Opt

using SkewLinearAlgebra
using Zygote
using LinearAlgebra

export Adam, step!, foo, SG, SGReg, AdamSpecial, Adam_multi, rg_update

foo() = println("foo")

# Struct containing all necessary info
mutable struct Adam{T<:Number}
    theta::Matrix{T} # Parameter array
    loss::Function                # Loss function
    m::SkewHermitian{T, Matrix{T}}     # First moment
    v::Hermitian{T, Matrix{T}}     # Second moment
    b1::Float64                   # Exp. decay first moment
    b2::Float64                   # Exp. decay second moment
    a::Float64                    # Step size
    eps::Float64                  # Epsilon for stability
    t::Int                  # Time step (iteration)
end

# Outer constructor
function Adam(theta::AbstractArray{T}, loss::Function) where T<:Number
    m = skewhermitian(zero(theta))
    v = Hermitian(zero(theta))
    b1 = 0.9
    b2 = 0.999
    a = 0.001
    eps = 1e-8
    t = 0
    Adam{T}(theta, loss, m, v, b1, b2, a, eps, t)
end

mutable struct Adam_multi{T<:Number, N}
    theta::NTuple{N, Matrix{T}}
    loss::Function                # Loss function
    m::NTuple{N, SkewHermitian{T, Matrix{T}}}     # First moment
    v::NTuple{N, Hermitian{T, Matrix{T}}}     # Second moment
    b1::Float64                   # Exp. decay first moment
    b2::Float64                   # Exp. decay second moment
    a::Float64                    # Step size
    eps::Float64                  # Epsilon for stability
    t::Int                  # Time step (iteration)
end

# Outer constructor
function Adam_multi(theta::NTuple{N, AbstractArray{T}}, loss::Function) where {N, T<:Number}
    m = ntuple(i -> skewhermitian(zero(theta[i])), Val(N))
    v = ntuple(i -> Hermitian(zero(theta[i])), Val(N))
    b1 = 0.9
    b2 = 0.999
    a = 0.001
    eps = 1e-8
    t = 0
    Adam_multi{T, N}(theta, loss, m, v, b1, b2, a, eps, t)
end

function rg_update(X::AbstractArray{T}, rg::AbstractArray{T}) where T<:Number
    return exp(-rg) * X
end

# function rg_update(X::AbstractArray{T}, rg::SkewHermitian{T, Matrix{T}}) where T<:Number
#     return exp(-rg) * X
# end

# Step function with optional keyword arguments for the data passed to grad()
function step!(opt::Adam{Float64})
    opt.t += 1
    gt′ = Zygote.gradient(opt.loss, opt.theta)[1]
    gt = skewhermitian((gt′ - gt′') / 2)
    opt.m = opt.b1 * opt.m + (1 - opt.b1) * gt
    opt.v = Hermitian(opt.b2 * opt.v + (1 - opt.b2) * (gt .^ 2))
    mhat = opt.m ./ (1 - opt.b1^opt.t)
    vhat = opt.v ./ (1 - opt.b2^opt.t)
    rg′ = opt.a .* (mhat ./ (sqrt.(vhat) .+ opt.eps))
    opt.theta = rg_update(opt.theta, rg′)
end


function step!(opt::Adam_multi{Float64})
    opt.t += 1
    gt′ = Zygote.gradient(opt.loss, opt.theta)[1]
    for i in 1:length(opt.theta)
        gt = skewhermitian((gt′[i] - gt′[i]') / 2)
        println(gt)
        opt.m[i][:] = opt.b1 * opt.m[i] + (1 - opt.b1) * gt
        opt.v[i][:] = Hermitian(opt.b2 * opt.v[i] + (1 - opt.b2) * (gt .^ 2))
        mhat = opt.m[i] ./ (1 - opt.b1^opt.t)
        vhat = opt.v[i] ./ (1 - opt.b2^opt.t)
        rg′ = opt.a .* (mhat ./ (sqrt.(vhat) .+ opt.eps))
        opt.theta[i][:] = rg_update(opt.theta[i], rg′)
    end
end

function step!(opt::Adam{Complex{Float64}})
    opt.t += 1
    gt′ = Zygote.gradient(opt.loss, opt.theta)[1]
    gt = skewhermitian((gt′ - gt′') / 2)
    opt.m = opt.b1 * opt.m + (1 - opt.b1) * gt
    vr = Hermitian(opt.b2 * real(opt.v) + (1 - opt.b2) * (real(gt) .^ 2))
    vi = Hermitian(opt.b2 * imag(opt.v) + (1 - opt.b2) * (imag(gt) .^ 2))
    opt.v = Hermitian(vr +  vi * 1im)
    mhat = opt.m ./ (1 - opt.b1^opt.t)
    mr = real(mhat)
    mi = imag(mhat)
    vhatr = vr ./ (1 - opt.b2^opt.t)
    vhati = vi ./ (1 - opt.b2^opt.t)
    rg′_r = opt.a .* (mr ./ (sqrt.(vhatr) .+ opt.eps))
    rg′_i = opt.a .* (mi ./ (sqrt.(vhati) .+ opt.eps))
    rg = rg′_r + 1im * rg′_i
    rg = skewhermitian(rg)
    opt.theta = rg_update(opt.theta, rg)
end

mutable struct AdamSpecial{T<:Number}
    theta::Matrix{T} # Parameter array
    loss::Function                # Loss function
    m::Matrix{T}     # First moment
    v::Matrix{T}     # Second moment
    b1::Float64                   # Exp. decay first moment
    b2::Float64                   # Exp. decay second moment
    a::Float64                    # Step size
    eps::Float64                  # Epsilon for stability
    t::Int                  # Time step (iteration)
end

function AdamSpecial(theta::AbstractArray{T}, loss::Function) where T<:Number
    m = zero(theta)
    v = zero(theta)
    b1 = 0.9
    b2 = 0.999
    a = 0.001
    eps = 1e-8
    t = 0
    AdamSpecial{T}(theta, loss, m, v, b1, b2, a, eps, t)
end

function step!(opt::AdamSpecial{Float64})
    opt.t += 1
    gt = Zygote.gradient(opt.loss, opt.theta)[1]
    opt.m = opt.b1 * opt.m + (1 - opt.b1) * gt
    opt.v = opt.b2 * opt.v + (1 - opt.b2) * (gt .^ 2)
    mhat = opt.m ./ (1 - opt.b1^opt.t)
    vhat = opt.v ./ (1 - opt.b2^opt.t)
    rg′ = opt.a .* (mhat ./ (sqrt.(vhat) .+ opt.eps))
    opt.theta = rg_update(opt.theta, rg′)
end

function step!(opt::AdamSpecial{Complex{Float64}})
    opt.t += 1
    gt = Zygote.gradient(opt.loss, opt.theta)[1]
    # gt = skewhermitian((gt′ - gt′') / 2)
    opt.m = opt.b1 * opt.m + (1 - opt.b1) * gt
    vr = opt.b2 * real(opt.v) + (1 - opt.b2) * (real(gt) .^ 2)
    vi = opt.b2 * imag(opt.v) + (1 - opt.b2) * (imag(gt) .^ 2)
    opt.v = vr +  vi * 1im
    mhat = opt.m ./ (1 - opt.b1^opt.t)
    mr = real(mhat)
    mi = imag(mhat)
    vhatr = vr ./ (1 - opt.b2^opt.t)
    vhati = vi ./ (1 - opt.b2^opt.t)
    rg′_r = opt.a .* (mr ./ (sqrt.(vhatr) .+ opt.eps))
    rg′_i = opt.a .* (mi ./ (sqrt.(vhati) .+ opt.eps))
    rg = rg′_r + 1im * rg′_i
    opt.theta = rg_update(opt.theta, rg)
end


mutable struct SG
    theta::Matrix{Float64} # Parameter array
    loss::Function         # Loss function
    a::Float64             # Step size
    t::Int                 # Time step (iteration)
end

# Outer constructor
function SG(theta::AbstractArray{Float64}, loss::Function; a::Float64=0.001)
    SG(theta, loss, a, 0)
end

function step!(opt::SG)
    opt.t += 1
    gt′ = Zygote.gradient(opt.loss, opt.theta)[1]
    gt = skewhermitian((gt′ - gt′') / 2) # Ensure the gradient is skew-Hermitian
    rg = opt.a * gt
    opt.theta = rg_update(opt.theta, rg)
end

mutable struct SGReg
    theta::Matrix{Float64} # Parameter array
    loss::Function         # Loss function
    reg::Function          # Regularization function
    a::Float64             # Step size
    t::Int                 # Time step (iteration)
    decay::Float64             # Decay rate
end

function SGReg(theta::AbstractArray{Float64}, loss::Function, reg::Function; a::Float64=0.001, decay::Float64=300.0)
    SGReg(theta, loss, reg, a, 0, decay)
end

function step!(opt::SGReg)
    opt.t += 1
    weight = exp(- opt.t / opt.decay)
    gt′ = Zygote.gradient(x -> (1 - weight) * opt.loss(x) + weight * opt.reg(x), opt.theta)[1]
    gt = skewhermitian((gt′ - gt′') / 2) # Ensure the gradient is skew-Hermitian
    rg = opt.a * gt
    opt.theta = rg_update(opt.theta, rg)
end

end

end # module MinEigLoss
