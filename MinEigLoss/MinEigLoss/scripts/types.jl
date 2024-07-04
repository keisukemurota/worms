function f(;t::Type{T}) where {T}
    T
end

g(::Type{T}) where {T} = T.var.ub

g(Array{<:Int})

eltype(Vector{Int})

f(t = Int)

Array{Int128} <: Array{<:Integer} 

Int128 <: Real

Integer |> subtypes

Int128 |> supertypes
Int |> supertypes   

isa(1.0, Integer)

Int128 <: Integer   

typejoin(Int, Int128)
typeintersect(Int, Int128)

isabstracttype(Int)
isconcretetype(Int)

typejoin(Int, Int128) 


f1(A::Array) = 1
f2(A::Array{Int}) = 2
f3(A::Array{T}) where {T<:Any} = 3
f4(A::Array{Any}) = 4

typeof(f3)


dump(Array)


T = TypeVar(:T,Integer)
UnionAll(T, Array{T})

supertypes(Int)

i :: Int32 = 1

Int64 <: typeof(Real(i))

Missing |> supertypes

[1, 2, 3] isa Vector{>:Int}
x = [1, missing, 3]
x isa Vector{>:Missing}
Union{Missing, Int} >: Missing

Union{Missing, Int} >: Int64
Missing |> supertypes

Vector{Missing} <: Vector{Any}

vv :: Vector{Union{Missing, Int}} = [1, missing, 3]
eltype(vv)

vv[1] |> typeof

