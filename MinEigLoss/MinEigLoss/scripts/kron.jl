using Pkg
Pkg.activate("./MinEigLoss/MinEigLoss")
using MinEigLoss
using SparseArrays
using Arpack
using Zygote
using Zygote: @adjoint
using LinearAlgebra
using EDKit
using BenchmarkTools

d = 3
u = RandomMatrix.Haar(d, Float64)

# function n_kron(u::AbstractMatrix{T}, ::Val{n}) where {T<:Real,n}
#     us = ntuple(i -> u, Val(n))
#     kron(us...)
# end

# n_kron(u) = n_kron(u, Val(5))

# loss(u) = norm(n_kron(u) - I)

# Zygote.gradient(loss, u)

# y, back = Zygote.pullback(n_kron, u)

# back(y)[1]

# randn(d^2, d^2) |> vec

# function d_n_kron(ȳ, u, ::Val{n}) where {n}
#     g = zero(u)
#     ukron = kron(ntuple(i -> u, Val(n - 1))...) * n
#     d = size(u, 1)
#     s = size(ukron, 1)
#     for i in 1:d
#         for j in 1:d
#             for k in 1:s
#                 for l in 1:s
#                     g[i, j] += ukron[k, l] * ȳ[(i-1)*s+k, (j-1)*s+l]
#                 end
#             end
#         end
#     end
#     g
# end

# using BenchmarkTools
# d = 3
# y = randn(d, d)
# y = n_kron(y, Val(5))

# d_n_kron(y, u, Val(5))
# back(y)[1]
# @benchmark d_n_kron($y, $u, Val(5))
# @benchmark back($y)[1]
# @code_warntype d_n_kron(y, u, Val(5))



using Arpack
# λ, v =  sparse(H) |> x -> eigs(x, which = :SR)

# function unitary_transform(H, u, ::Val{n}) where {n}
#     ukron = n_kron(u, Val(n))
#     ukron * H * ukron'
# end

# Harray = H |> Array

# @benchmark unitary_transform($Harray, $u, Val(8))

# @benchmark sparse($Harray)

# @benchmark sparse($h)



function column_mod!(target::AbstractVecOrMat, M::SparseMatrixCSC, I::Vector{Int}, b::AbstractBasis, coeff=1)
    rows, vals = rowvals(M), nonzeros(M)
    j = index(b.dgt, I, base=b.B)
    change = false
    for i in nzrange(M, j)
        row, val = rows[i], vals[i]
        change!(b.dgt, I, row, base=b.B)
        C, pos = index(b)
        isa(target, AbstractVector) ? (target[pos] += coeff * C * val) : (target[pos, :] .+= (C * val) .* coeff)
        change = true
    end
    change && change!(b.dgt, I, j, base=b.B) 
    nothing
end
#---------------------------------------------------------------------------------------------------
function colmn_mod!(target::AbstractVecOrMat, opt::EDKit.Operator, j::Integer, h :: SparseMatrixCSC, coeff=1)
    b, M, I = opt.B, opt.M, opt.I
    r = change!(b, j)
    C = isone(r) ? coeff : coeff / r
    hs = sparse(h)
    for i = 1:length(M)
        column_mod!(target, hs, I[i], b, C)
    end
end


function trans_inv_ham(h :: AbstractArray, opt :: EDKit.Operator)
    h = sparse(h)
    H = spzeros(eltype(opt), size(opt)...)
    if size(opt.M, 1) > 0 && size(opt.M, 2) > 0
        for j = 1:size(opt.B, 2)
            colmn_mod!(view(H, :, j), opt, j, h)
        end
    end
    H
end

function adjoint_sys_ham(h :: AbstractArray, sys_grad)
    loc_s = size(h, 1)
    unit_s = sqrt(size(h, 1))
    sys_s = size(sys_grad, 1)
    n_site = log(unit_s, sys_s)
    rem_s = div(sys_s, loc_s)
    grad = zero(h)
    sys_grad = sys_grad .* n_site
    dh = zero(h)
    for k = 1:loc_s
        for l = 1:loc_s
            for i = 1:rem_s
                grad[k, l] += sys_grad[rem_s * (k - 1) + i, rem_s * (l - 1) + i]
            end
        end
    end
    grad
end

@adjoint trans_inv_ham(h, opt) = trans_inv_ham(h, opt), c -> (adjoint_sys_ham(h, c), nothing)

function loss_func(h, opt)
    H = sparse(opt)
    # e0 = -eigs(H, which = :SR)[1][1]
    function wrapper(u)
        u = MinEigLoss.unitary(u)
        U = kron(u, u)
        h̄ = U * h * U'
        h̄ = (h̄ + h̄') ./ 2
        H = trans_inv_ham(h̄, opt)
        H = abs.(H)
        ϕ = zeros(size(H, 1))
        Zygote.ignore() do 
            ϕ .= eigs(H, nev = 1, which = :LR)[2][:, 1]
        end
        ϕ'*H*ϕ
    end
    wrapper
end

d = 3; # spin degree of freedom
χ = 2; # number of bond dimensions
h = FF.ff(d, χ, dim = Val(1), seed = 1);
L = 6;
h = h - eigmax(h)I;
opt = trans_inv_operator(h, 1:2, L) ;

loss_u = loss_func(h, opt);
mle_sys_unitary = MinEigLoss.mle_sys_unitary(opt |> Array, d);

gu = Zygote.gradient(loss_u, u)[1]
gu2 = Zygote.gradient(mle_sys_unitary, u)[1]
u2 = MinEigLoss.Opt.rg_update(u, gu * 0.01)
mle_sys_unitary(u2) - mle_sys_unitary(u)
loss_u(u2) - loss_u(u)



@benchmark loss_u($u)
@benchmark mle_sys_unitary($u)
@benchmark Zygote.gradient($loss_u, $u)[1]
@benchmark Zygote.gradient(mle_sys_unitary, u)[1]


Zygote.gradient(2, 3) do a, b
    Zygote.@showgrad(a) * b
    # a * Zygote.@showgrad(b)
end

# U = kron(u, u, u, u)
# H = opt |> Array
# U' * H * U - loss_u(u)

# loss_u(u)|>Array





# h = trans_inv_ham(h, opt)
# gh

# sign.(h) - sign.(gh)

# H = trans_inv_ham(h, opt)

# H_array = H |> Array

# H_mat = reshape(H_array, ntuple(i -> d, Val(2*L))...)





# function trans_inv_1d(h :: AbstractArray, L :: Int)
#     opt = trans_inv_operator(h, 1:2, L)
#     M = opt.M
#     I = opt.I

#     @adjoint create_sys_ham(h) = create_sys_ham(h), c -> (adjoin_sym_ham(h, c), )

#     function create_sys_ham(h :: AbstractArray)
#         h = sparse(h)
#         H = spzeros(eltype(opt), size(opt)...)
#         if size(M, 1) > 0 && size(M, 2) > 0
#             for j = 1:size(opt.B, 2)
#                 colmn_mod!(view(H, :, j), opt, j, h)
#             end
#         end
#         H
#     end
# end

system_ham = trans_inv_1d(h, L)


@code_warntype system_ham(h)


@benchmark unitary_transform($opt, $u)
sparseH
eigs(sparseH)
sparse(opt) |> typeof
sparseH |> typeof
eigs(sparse(opt))
λ, v = eigs(sparseH)

λ

sparse(opt)

sparseH |> sum
nonzeros(sparseH - sparseH')




λ, v = eigs(M, which = :SR, nev = 1)

@benchmark eigs($M, which = :SR, nev = 1)


v = ones(size(M, 1))
v ./= norm(v)
@benchmark $M * $v

h_ = h |> u -> view(u, :, 1)
h[1] = 3
# @code_warntype trans_inv_operator(h, 1:2, L)
# @code_typed trans_inv_operator(h, 1:2, L)

# 1:2 |> Array


# base = EDKit.find_base(size(h, 1), length(1:2))
# B = TensorBasis(L=L, base=base)
# # @code_warntype trans_inv_operator(h, 1:2, B)
# # @which trans_inv_operator(h, 1:2, B)



# ind = 1:2
# L = length(B.dgt)
# smat = sparse(h)
# mats = fill(smat, L)
# inds = [mod.(ind .+ i, L) .+ 1 for i = -1:L-2]

# function operator_from_mat(mats::AbstractVector{<:AbstractMatrix}, inds::AbstractVector{<:AbstractVector}, B::AbstractBasis)
#     num = length(mats)
#     @assert num == length(inds) "Numbers mismatch: $num matrices and $(length(inds)) indices."
#     dtype = promote_type(eltype.(mats)...)
#     M = Vector{SparseMatrixCSC{dtype, Int64}}(undef, num)
#     I = Vector{Vector{Int64}}(undef, num)
#     N = 0
#     for i = 1:num
#         iszero(mats[i]) && continue
#         ind = inds[i]
#         pos = findfirst(x -> isequal(x, ind), view(I, 1:N))
#         if isnothing(pos)
#             N += 1
#             I[N] = ind
#             M[N] = sparse(mats[i])
#         else
#             M[pos] += mats[i]
#         end
#     end
#     deleteat!(M, N+1:num)
#     deleteat!(I, N+1:num)
#     EDKit.Operator(M, I, B)
# end

# H = operator_from_mat(mats, inds, B)

# @benchmark sparse(H)






