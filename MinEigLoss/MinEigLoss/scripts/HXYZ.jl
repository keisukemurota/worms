using Pkg
Pkg.activate("./MinEigLoss/MinEigLoss")
using MinEigLoss
using Revise
using LinearAlgebra
using EDKit
using Arpack
using Test
using ProgressBars
using Printf
using Random
using JLD2
using PyCall
using Zygote

np = pyimport("numpy")


function mle_sys_unitary_free(H:: AbstractMatrix, s_loc :: Int)
    s_sys = size(H, 1)
    n_rep :: Int = log(s_loc, s_sys) |> round |> Int
    e0 = eigmin(H)
    function kron_w(ws) :: AbstractMatrix
        @assert n_rep > 1
        if n_rep != length(ws)
            throw(ArgumentError("length of ws must be $n_rep"))
        end
        W = kron(ws[1], ws[2])
        
        for i in 3:n_rep
            W = kron(W, ws[i])
        end
        return W
    end

    function wrapper_special(ws)
        ws_unitary = ntuple(i -> MinEigLoss.unitary(ws[i]), length(ws))
        W = kron_w(ws_unitary)
        H̄ = (W * H * W') .|> abs
        return (eigvals(H̄) |> real |> maximum) + e0
    end

    return wrapper_special
end


import Base: inv

function inv(v::Tuple{Vararg{AbstractMatrix}})
    return map(inv, v)
end

function mle_sys_sim_free(H:: AbstractMatrix, s_loc :: Int)
    s_sys = size(H, 1)
    n_rep :: Int = log(s_loc, s_sys) |> round |> Int
    e0 = eigmin(H) |> abs
    function kron_v(vs) :: AbstractMatrix
        @assert n_rep > 1
        if n_rep != length(vs)
            throw(ArgumentError("length of ws must be $n_rep"))
        end
        V = kron(vs[1], vs[2])
        
        for i in 3:n_rep
            V = kron(V, vs[i])
        end
        return V
    end

    function wrapper_special(ws)
        ws_unitary = ntuple(i -> MinEigLoss.unitary(ws[i]), length(ws))
        W = kron_w(ws_unitary)
        W_inv = kron_v(inv(ws))
        H̄ = (W * H * W_inv) .|> abs
        return (eigvals(H̄) |> real |> maximum) + e0
    end

    return wrapper_special
end


function random_theta(d, type = :orthogonal)
    if type == :orthogonal
        return RandomMatrix.Haar(d, Float64)
    elseif type == :unitary
        return RandomMatrix.Haar(d, ComplexF64)
    else
        throw(ArgumentError("Invalid type: $type"))
    end
end


function heisenberg(;Js :: Vector{<:Real}, hx :: Float64 = 0.0, dim :: Int = 1)

    @assert length(Js) == 3

    h = (spin("XX") * Js[1] + spin("YY") * Js[2] + spin("ZZ") * Js[3]) ./ 4
    h = h + (hx / 2dim) .* (spin("XI") + spin("IX")) ./2
    return h
end

function lattice(; h :: AbstractMatrix{<:Real}, Ls :: Vector{<:Int}, BC :: Symbol = :PBC)

    @assert length(Ls) > 0

    if length(Ls) == 1
        bonds = BC == :PBC ? [[i, mod(i, Ls[1]) + 1] for i in 1:Ls[1]] : [[i, i+1] for i in 1:Ls[1]-1]
        println(bonds)
        op = EDKit.operator(fill(h, length(bonds)), bonds, Ls[1])
        return op, bonds
    else
        throw(ArgumentError("Invalid length of Ls: $(length(Ls))"))
    end

    return h
end

h = heisenberg(Js = [1.0, 1.0, 0.0], hx = 0.0, dim = 1);
offset = eigmax(h |> Array)
h = h - offset * I
h⁺ = -(h .|> abs)


Ls = [7] # if lattice size is odd, it has negative sign problem
op, bonds = lattice(h = h, Ls = Ls, BC = :PBC)
op⁺, bonds⁺ = lattice(h = h⁺, Ls = Ls, BC = :PBC)
H = op |> EDKit.sparse
H⁺ = op⁺ |> EDKit.sparse

λ, ϕ = eigs(H + offset * I * length(bonds), which=:SR, nev = 3)
λ⁺, ϕ⁺ = eigs(H⁺ + offset * I * length(bonds⁺), which=:SR, nev = 3)

@test λ ≈ λ⁺



loss = mle_sys_unitary_free(H |> Array |> Symmetric, 2)

loss(ntuple(i -> I(2), Val(prod(Ls))))

begin epsilon = 0.05; ws = ntuple(i -> random_theta(2, :orthogonal), Val(prod(Ls)))
    for n in 1:100
        g = Zygote.gradient(loss, ws)[1]

        for i in 1:prod(Ls)
            ws[i][:] = Opt.rg_update(ws[i], g[i] * epsilon)
        end
    end
    loss(ws)
end

# A |> eigvals

# A .|> abs |> eigvals

# V |> eigvals

# V |> eigvals

# (B = randn(4, 4)) |> eigvals

# B |> eigvals

# (B |> eigen).vectors

# exp(B)


function spec_test(w)
    global h
    h′ = w * h * inv(w)
    return (h′ |> diag |> minimum)
end

h = randn(4, 4)
h = h + h'

h = -(h - eigmax(h) * I)
# h = -(h - 100 * I)



spec_test(I(4))

begin w = randn(4, 4); epsilon = 0.01
    for i in 1:100
        g = Zygote.gradient(spec_test, w)[1]
        w = w - g * epsilon
        w = w ./ ((det(w) |> abs) ^ (1/4)) 
        @show spec_test(w)
    end
end

h′ = w * h * inv(w)

w 

h′

(w ./ ( (det(w) |> abs) ^ (1/4))) |> det








eigvals(h)
h

h′ |> eigvals

h′


spec_test(h1 |> Array, I(4))


u = ws[1] |> MinEigLoss.unitary

u

kron(u, u)





# MinEigLoss.rg_update(I(2), g[1])

# adam_muti = Opt.Adam_multi(ws, loss)

# Opt.step!(adam_muti)

function unitary_optimizer(d, optimizing_loss, tracking_loss; epochs, runs, params, type = :orthogonal)
    best_thetas = []
    best_losses_opt = []
    best_losses_track = []
    function random_theta()
        if type == :orthogonal
            return RandomMatrix.Haar(d, Float64)
        elseif type == :unitary
            return RandomMatrix.Haar(d, ComplexF64)
        else
            throw(ArgumentError("Invalid type: $type"))
        end
    end
    
    for i in 1:runs
        theta = random_theta()
        adam = Opt.Adam(theta, optimizing_loss)
        for (key, value) in params
            setfield!(adam, key, value)
        end
        
        run_best_theta = nothing
        run_best_loss_opt = Inf
        run_best_loss_track = Inf
        
        for j in 1:epochs
            Opt.step!(adam)
            current_loss_opt = optimizing_loss(adam.theta)
            current_loss_track = tracking_loss(adam.theta)
            
            if current_loss_opt < run_best_loss_opt
                run_best_loss_opt = current_loss_opt
                run_best_theta = copy(adam.theta)
                run_best_loss_track = current_loss_track
            end
        end
        
        push!(best_thetas, run_best_theta)
        push!(best_losses_opt, run_best_loss_opt)
        push!(best_losses_track, run_best_loss_track)
    end
    
    return (thetas = best_thetas, losses_opt = best_losses_opt, losses_track = best_losses_track)
end


function main(;Js :: Vector{<:Real}, hx :: Float64 = 0.0, dim :: Int = 1)
    d = 2
    h = heisenberg(Js = Js, hx = hx, dim = dim) |> Array

    offset = eigmax(h)
    h = h - offset * I

    mle_unitary = MinEigLoss.mle_unitary(h)

    # println("local hamiltonian : ", h)
    # println("mle_unitary : ", mle_unitary(I(2)))

    # @show h;

    uni_loc = unitary_optimizer(d, mle_unitary, mle_unitary, epochs = 300, runs = 50, params = Dict(:a => 0.02, :b1 => 0.95), type = :unitary);
    # println(uni_loc.losses_opt |> minimum)
    # println(uni_loc.losses_track |> minimum)

    orth_loc = unitary_optimizer(d, mle_unitary, mle_unitary, epochs = 300, runs = 50, params = Dict(:a => 0.02, :b1 => 0.95), type = :orthogonal);
    # println(orth_loc.losses_opt |> minimum)
    # println(orth_loc.losses_track |> minimum)

    return (uni_loc = uni_loc, orth_loc = orth_loc, mle_unitary = mle_unitary, h = h)
end


begin 
    h1 = heisenberg(Js = [1, 1, 0], hx = 0.0, dim = 1);
    h1 = h1 - eigmax(h1 |> Array) * I;
    H1, bonds1 = lattice(h = h1, Ls = [7], BC = :PBC);
    H1 = H1 |> Array
    # h2 = heisenberg(Js = [1.0, +1.3, 1.0], hx = 0.0, dim = 1);
    # h2 = h2 - eigmax(h2 |> Array) * I;
    # H2, bonds2 = lattice(h = h2, Ls = [3], BC = :PBC);
end


H1

mle_sys_unitary = MinEigLoss.mle_sys_unitary(H1, 2)
loss_unitary = MinEigLoss.mle_unitary(h1 |> Array)
loss_unitary(I(2))

uni_sys = unitary_optimizer(2, mle_sys_unitary, loss_unitary, epochs = 300, runs = 50, params = Dict(:a => 0.02, :b1 => 0.95), type = :orthogonal);
uni_sys.losses_opt |> minimum
mle_sys_unitary(I(2))

uni_sys.losses_opt |> argmin

uni_sys.thetas[argmin(uni_sys.losses_opt)]

U = kron(u, u, u)

MinEigLoss.mle_sys_unitary(U * H1 * U', 2)(I(2))







H1 |> Array
loss1 = MinEigLoss.mle_unitary(h1)
loss1(I(2))
loss1(random_theta(2))


h2  |> eigvals
(-abs.(h2)) |> eigvals


h2




res = main(Js = [1.3, 2, 1.0], hx = 0.0, dim = 1);
res.h
res = main(Js = [1.0, +1.3, 1.0], hx = 0.0, dim = 1);
res.h 
res.h |> eigvals
-(res.h .|> abs) |> eigvals
println(res.uni_loc.losses_opt |> minimum)
println(res.orth_loc.losses_opt |> minimum)
println(res.mle_unitary(I(2)))

u = exp(im * π / 2 .* spin("Z") |> Array)

u * spin("Y") * u' + spin("Y")

u'

@show res.h



# Print the working directory
h = -np.load("python/rmsKit/array/torch/HXYZ1D_loc/Jx_1.3_Jy_2_Jz_1_hx_0_hz_0/1_mel/H/0.npy")
u = -np.load("python/rmsKit/array/torch/HXYZ1D_loc/Jx_1.3_Jy_2_Jz_1_hx_0_hz_0/1_mel/Adam/lr_0.01_epoch_100/loss_0.0000024/u/0.npy")

h = h - eigmax(h) * I


h

eigvals(h)
eigvals(-(h .|> abs))

U = kron(u, u)



loss = MinEigLoss.mle_unitary(h)

loss(I(2))
loss(u)
