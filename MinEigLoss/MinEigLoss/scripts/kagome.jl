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

function random_theta(d, type = :orthogonal)
    if type == :orthogonal
        return RandomMatrix.Haar(d, Float64)
    elseif type == :unitary
        return RandomMatrix.Haar(d, ComplexF64)
    else
        throw(ArgumentError("Invalid type: $type"))
    end
end


function unitary_optimizer(d, optimizing_loss, tracking_loss; epochs, runs, params :: Dict, type = :orthogonal)
    best_thetas = []
    best_losses_opt = []
    best_losses_track = []
    
    iter = ProgressBar(1:runs)
    for i in iter
        theta = random_theta(d, type)
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
        set_description(iter, string(@sprintf("Loss: %.2f", best_losses_opt |> minimum)))
    end
    
    return (thetas = best_thetas, losses_opt = best_losses_opt, losses_track = best_losses_track)
end

function kagome_heisenberg(;Js :: Vector{<:Real}, hx :: Float64 = 0.0)

    # Js = [Jx, Jy, Jz]
    @assert length(Js) == 3

    h = spin("XX") * Js[1] + spin("YY") * Js[2] + spin("ZZ") * Js[3]
    h = h ./ 4

    # hx = hx * spin("X")

    bonds = [[1, 2], [2, 3], [3, 1], [4, 5], [5, 6], [6, 4]]
    mats = fill(h, 6)
    h̄ = EDKit.operator(mats, bonds, 6) |> Array

    hs = fill(h̄, 3)

    hs[1] += EDKit.operator([h], [[3, 4]], 6) |> Array  
    hs[2] += EDKit.operator([h], [[2, 4]], 6) |> Array
    hs[3] += EDKit.operator([h], [[3, 5]], 6) |> Array

    return hs
end

Jx = -2:0.2:2
Jy = -2:0.2:2

res = Dict{Tuple{Float64, Float64}, Any}()

for jx in Jx, jy in Jy

    if jy <= -jx
        continue
    end
    J1 = [jx, jy, 1]
    hs = kagome_heisenberg(Js = J1, hx = 0.0)
    hs = [h - eigmax(h)I for h in hs]

    mle_unitary_losses = [MinEigLoss.mle_unitary(h) for h in hs]
    mle_unitary_loss(u) = sum([loss(u) for loss in mle_unitary_losses])

    println("jx: $jx, jy: $jy")
    res_o = unitary_optimizer(8, mle_unitary_loss, x -> nothing, epochs = 400, runs = 20, params = Dict(:a => 0.02, :b1 => 0.95), type = :orthogonal);

    minidx = findmin(res_o.losses_opt)[2]
    res[(jx, jy)] = (loss = res_o.losses_opt[minidx], theta = res_o.thetas[minidx], h = hs)
end

res = load_object(joinpath(dirname(@__FILE__),"pickles", "KH_unitary.jld2"))
# load


z = Float64[]
function color(jx, jy)
    if jy <= -jx
        return 0.0
    end
    hs = res[(jx, jy)].h
    mle_unitary_losses = [MinEigLoss.mle_unitary(h) for h in hs]
    mle_unitary_loss(u) = sum([loss(u) for loss in mle_unitary_losses])

    # e0 = mle_unitary_loss(I(8))
    e1 = res[(jx, jy)].loss / mle_unitary_loss(I(8))
    return e1
end


# using Plots
heatmap(
    Jx, 
    Jy, 
    [color(jx, jy) for jy in Jy, jx in Jx], 
    size = (700, 700),
    aspect_ratio = :equal,
    xlims = (-2, 2),
    ylims = (-2, 2),
    xlabel = "Jx",
    ylabel = "Jy",
    title = "L1 adaptive loss (unitary)",
    colorbar = true,
    c = :Reds,
    # margin = 0Plots.mm,
    rightmargin = 2Plots.mm,
    titlefontsize = 20,
    guidefontsize = 12,
    tickfontsize = 12,
    legendfontsize = 12,
)



# mle_unitary_loss(I(8))
# mle_unitary_loss(random_theta(8, :unitary))

# using BenchmarkTools
# @benchmark mle_unitary_loss($I(8))


# res_o

# begin type=:orthogonal;  optimizing_loss = mle_unitary_loss; epochs = 700; params = Dict(:a => 0.01, :b1 => 0.95); d = 8
#     losses = []
#     p = plot(title="Optimization Progress", xlabel="Iteration", ylabel="Loss", legend=false)
#     theta = random_theta(d, type)
#     adam = Opt.Adam(theta, optimizing_loss)
#     seed = rand(1:1000000)
#     println(seed)
#     Random.seed!(seed)
#     for (key, value) in params
#         setfield!(adam, key, value)
#     end
#     iter = ProgressBar(1:epochs)
#     for j in iter
#         Opt.step!(adam)
#         current_loss_opt = optimizing_loss(adam.theta)
#         set_description(iter, string(@sprintf("Loss: %.2f", current_loss_opt)))
#         push!(losses, current_loss_opt)
#     end
#     plot!(p, losses, ylim=(0, maximum(losses)), show=true)
#     display(p)
# end


# res_o.losses_opt |> minimum
# res_o.losses_opt |> minimum
