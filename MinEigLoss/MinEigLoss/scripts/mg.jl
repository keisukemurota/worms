using Pkg
Pkg.activate("./MinEigLoss/MinEigLoss")
using MinEigLoss
using Revise
using LinearAlgebra
using EDKit
using Arpack
using Test
using SparseArrays


function mg_model(Js :: Vector{<:Real}, L_half :: Int)
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
    H = Array(op)
    Es = eigvals(H)
    for (i, β) in enumerate(βs)
        boltz = exp.(-β .* Es)
        Z = sum(boltz)
        E = sum(boltz .* Es) / Z
        C = β^2 * (sum(boltz .* Es .^ 2) / Z - E^2)
        push!(res, (β=β, E=E/L, C=C/L))
    end
    return res
end

L_half = 5
L = L_half * 2
e_mg = -L * (3/4.0) 

J_mg = [1, 2, 2]

op = mg_model(J_mg, L_half);
obs = observable(op, [1.0, 2.0, 4.0]);

v = map(obs) do x
    (β=x.β, E=x.E * 2, C=x.C * 2)
end;

Hs = op |> EDKit.sparse;
λ_mg, ϕ = eigs(Hs, which=:SR, nev = 3)
@test λ_mg[1] ≈ e_mg



op = mg_model([1, 2, 1], L_half);
Hs = op |> EDKit.sparse;
λ1, _ = eigs(Hs, which=:SR, nev = 3)
op = mg_model([1, 1, 2], L_half);
Hs = op |> EDKit.sparse;
λ2, _ = eigs(Hs, which=:SR, nev = 3)

@test λ1 ≈ λ2




function mg_local(Js :: Vector{<:Real}, ::Val{:lt1})
    # J1, J2, J3 model
    h = spin("XX") + spin("YY") + spin("ZZ")
    h = h ./ 4

    @assert length(Js) == 3

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


h = mg_local([1, 2, 2], Val(:lt1))

h |> EDKit.Array

H = EDKit.trans_inv_operator(h, 1:2, L_half) |> EDKit.sparse

λ_dimer, ϕ = eigs(H, which=:SR, nev = 3)

@test λ_dimer ≈ λ_mg


u = [
    1 1 0 0;
    0 0 1 1;
    0 0 1 -1;
    -1 1 0 0;
] ./ sqrt(2)

U = kron(u, u)

begin
    h = -mg_local([1, 2, 2], Val(:lt1)) |> Array
    h -= I(size(h,1))*eigmin(h)
    h  = U' * h * U |> Symmetric |> sparse
    h⁺ = abs.(h)


    @test diag(h) ≈ diag(h⁺)

    bonds = [[i, mod(i, L_half) + 1] for i in 1:L_half-1]
    op = EDKit.operator(fill(h, L_half-1), bonds, L_half)
    op⁺ = EDKit.operator(fill(h⁺, L_half-1), bonds, L_half) 
    λ = eigmax(Array(op))
    λ⁺ = eigmax(Array(op⁺))

    @test λ ≈ λ⁺
end
# op = EDKit.trans_inv_operator(h, 1:2, L_half)
# op⁺ = EDKit.trans_inv_operator(h⁺, 1:2, L_half)

