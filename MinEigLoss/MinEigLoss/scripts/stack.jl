using LinearAlgebra
using SparseArrays
using IterativeSolvers



N = 6000
A = sprandn(N, N, N / N^2)
A = A * A'
b = randn(N)
x = zero(b)
A \ b

using Arpack

eigs(A, nev=1, maxiter=10, which=:LR)[5]

IterativeSolvers.cg(A, b)

@code_typed mul!(x, A, b)

randn(10)

Base.summarysize(Float64(1.0))

