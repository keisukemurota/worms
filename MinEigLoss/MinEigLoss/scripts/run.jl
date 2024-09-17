using Pkg
Pkg.activate("./MinEigLoss/MinEigLoss")
using Random
using MinEigLoss
using Revise
using LinearAlgebra
using EDKit
using Plots
using Zygote
using ITensors
using JLD2

RandomMatrix.Haar(3, ComplexF64) |> det

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

RandomMatrix.Similarity(3, ComplexF64)

function special_optimizer(d, optimizing_loss, tracking_loss; epochs, runs, params, type = :orthogonal)
    best_thetas = []
    best_losses_opt = []
    best_losses_track = []
    function random_theta()
        if type == :orthogonal
            w = RandomMatrix.Similarity(d, Float64)
        elseif type == :unitary
            w = RandomMatrix.Similarity(d, ComplexF64)
        else
            throw(ArgumentError("Invalid type: $type"))
        end
        return w
    end
        
    for _ in 1:runs
        theta = random_theta()
        adam = Opt.AdamSpecial(theta, optimizing_loss)
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



function heisenberg(;Js :: Vector{<:Real}, hx :: Float64 = 0.0, dim :: Int = 1, D :: Int = 2)

    @assert length(Js) == 4

    h = spin((Js[1], "xx"), (Js[2], "yy"), (Js[3], "zz"), (Js[4], "zx"), D=D)
    h = h + (hx / 2dim) .* (spin((1, "x1"), D = D) + spin((1, "1x"), D = D)) ./2
    return h
end

heisenberg(Js = [1, 1, 1, 1], hx = 0.0, dim = 1, D = 3)


function mg_local(Js :: Vector{<:Real})
    @assert length(Js) == 3
    h = spin("XX") + spin("YY") + spin("ZZ")
    h = h ./ 4
    hs = map(Js) do J
        J .* h
    end
    bonds = [
        [1, 2], [3, 4],
        [2, 3],
        [1, 3], [2, 4]
    ]
    mats = [
        hs[2] ./ 2, hs[2] ./ 2,
        hs[3],
        hs[1], hs[1]
    ]
    EDKit.operator(mats, bonds, 4) |> EDKit.sparse
end

function main(seed ; Js = [], hx = [])
    d = 3; # spin degree of freedom
    χ = 2; # number of bond dimensions
    h = FF.ff(d, χ, dim = Val(1), seed = seed);
    # Generate random coupling constants
    # Random.seed!(seed);
    # if Js == []
    #     Js = rand(Float64, 3) .* 2 .- 1
    #     Js .*= 2
    #     hx = rand(Float64) *2 - 1
    # end
    print("Js : ", Js, "hx : ", hx)
    # h = heisenberg(Js = Js, hx = hx, dim = 1, D = d) |> Array
    # h = mg_local(Js) |> Array
    L = 3;
    # h = h - eigmax(h)I
    # h = h - 100I
    H = trans_inv_operator(h, 1:2, L) |> Array;

    mle_unitary = MinEigLoss.mle_unitary(h)
    mle_sys_unitary = MinEigLoss.mle_sys_unitary(H, d)
    mle_sys_special = MinEigLoss.mle_sys_special(H, d)
    l1_sys_special = MinEigLoss.l1_sys_special(H, d)
    mle_special = MinEigLoss.mle_special(h)

    if mle_sys_unitary(I(d)) < 1e-4
        return (spec_sys_orth = [], spec_sys_uni = [], uni_sys = [], orth_sys = [], h = h, H = H)
    end

    println("initial loss unitary: ", mle_sys_unitary(I(d)))
    println("initial loss special: ", mle_sys_special(I(d)))

    # uni_sys = []
    # uni_sys = unitary_optimizer(d, mle_sys_unitary, mle_unitary, epochs = 200, runs = 50, 
                                            # params = Dict(:a => 0.02, :b1 => 0.95), type = :unitary);
    # println(uni_sys.losses_opt |> minimum)

    orth_sys = unitary_optimizer(d, mle_sys_unitary, mle_unitary, epochs = 200, runs = 50, 
                                            params = Dict(:a => 0.03, :b1 => 0.95), type = :unitary);
    println(orth_sys.losses_opt |> minimum)
    # println(orth_sys.losses_track |> minimum)

    spec_sys_orth = special_optimizer(d, mle_sys_special, mle_special, epochs = 200, runs = 50, 
    params = Dict(:a => 0.03, :b1 => 0.95), type = :unitary);
    println(spec_sys_orth.losses_opt |> minimum)

    spec_sys_uni = []

    return (spec_sys_orth = spec_sys_orth, spec_sys_uni = spec_sys_uni, uni_sys = uni_sys, orth_sys = orth_sys, h = h, H = H)
    # return (spec_sys = spec_sys, h = h, H = H)
end


# h = heisenberg(Js = [1, 1, 1], hx = 0.0, dim = 1, D = 3) |> Array



res = Dict()
for i in 1:100
    seed = rand(10 ^ 5 : 2 * 10^6)  # Generate random seed between 1000 and 9999
    result = main(seed)
    res[seed] = result
    println("seed: $seed is done")
    println("--" ^ 20)
end

# r = main(320903);

# r2 = main(320903, Js = [-1, 0.5, 1.5, 0], hx = 0.0);

Jy = 0:0.2:4
Jz = 0:0.2:4
res = Dict()
for jy in Jy
    for jz in Jz
        res[jy, jz] = main(320903, Js = [1, jy, jz], hx = 0.0)
        println("jy: $jy, jz: $jz is done")
    end
end


save_object("./spec_vs_uni_complex.jld2", res)



heisenberg(Js = [-1.4043701602807586, 0.37542649015580265, 1.717261834125226, 0], hx = 0.729469172122322, dim = 1, D = 3) |> Array

res

for (key, value) in res
    if value.orth_sys != []
        # println(key)
        if value.orth_sys.losses_opt |> minimum < value.spec_sys.losses_opt |> minimum
            println("orth_sys is better")
        else
            println("spec_sys is better, ", key)
            println("orth_sys: ",value.orth_sys.losses_opt |> minimum)
            println("spec_sys: ",value.spec_sys.losses_opt |> minimum)
        end
    end
end
# save_object(joinpath(dirname(@__FILE__),"pickles", "res_300.jld2"), res)


# (res |> values |> collect)[1].h

d = 3; # spin degree of freedom
χ = 2; # number of bond dimensions
h = FF.ff(d, χ, dim = Val(1), seed = 456773);
L = 6;
h = h - eigmax(h)I
H = trans_inv_operator(h, 1:2, L) |> Array;

mle_unitary = MinEigLoss.mle_unitary(h)
mle_sys_unitary = MinEigLoss.mle_sys_unitary(H, d)
mle_sys_special = MinEigLoss.mle_sys_special(H, d)
l1_sys_special = MinEigLoss.l1_sys_special(H, d)
mle_special = MinEigLoss.mle_special(h)


uni_loc = unitary_optimizer(d, mle_unitary, mle_sys_unitary, epochs = 100, runs = 1, params = Dict(:a => 0.02, :b1 => 0.95), type = :unitary);

