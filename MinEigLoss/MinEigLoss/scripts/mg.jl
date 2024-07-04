using Pkg
Pkg.activate("./MinEigLoss/MinEigLoss")
using MinEigLoss
using Revise
using LinearAlgebra
using EDKit
using Arpack
using Test


function mg_model(Js :: Vector{Float64}, L_half :: Int)
    L = L_half * 2
    # J1, J2, J3 model
    h = spin("XX") + spin("YY") + spin("ZZ")
    h = h ./ 4

    @assert length(Js) == 3

    hs = map(Js) do J
        J .* h
    end

    bonds1 = [[i, mod(i + 1, L) + 1] for i in 1:L];
    bonds2 = [[2*i, mod(2*i, L) + 1] for i in 1:L_half];
    bonds3 = [[2*i - 1, mod(2*i-1, L) + 1] for i in 1:L_half];

    mats = [
        fill(hs[1], L); 
        fill(hs[2], L_half);
        fill(hs[3], L_half);
    ];

    bonds = [
        bonds1;
        bonds2;
        bonds3;
    ];

    EDKit.operator(mats, bonds, L);
end

restype = NamedTuple{(:β, :E, :C), Tuple{Float64, Float64, Float64}}

function observable(op, βs :: Vector{Float64})
    L = size(op.B.dgt)[1]
    res :: Vector{restype} = []
    for (i, β) in enumerate(βs)
        H = Array(op)
        Es = eigvals(H)
        boltz = exp.(-β .* Es)
        Z = sum(boltz)
        E = sum(boltz .* Es) / Z
        C = β^2 * (sum(boltz .* Es .^ 2) / Z - E^2)
        push!(res, (β=β, E=E/L, C=C/L))
    end
    return res
end

op = mg_model([1, 1.8, 2.4], 5);

obs = observable(op, [1.0, 2.0, 4.0]);

map(obs) do x
    (β=x.β, E=x.E * 2, C=x.C * 2)
end

Hs = op |> EDKit.sparse;

e_mg = -L * (3/4.0) * 1/2.0

λ, ϕ = eigs(Hs, which=:SR, nev = 1)

@test λ[1] ≈ e_mg



H = Array(op)
Es = eigvals(H)

beta = [1.0, 2.0]

Z = exp.(-beta .* reshape(Es, 1, :)) |> (e -> sum(e, dims=2))

sum(beta)

reshape(beta, 1, :) .* Es
