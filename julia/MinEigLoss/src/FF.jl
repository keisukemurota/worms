module FF
    export ff, hello
    hello() = println("Hello, World!")
    using Random
    using ITensors
    function ff(d :: Int64, χ :: Int64; dim :: Val{1}, seed :: Union{Int64, Nothing} = nothing)
        if !isnothing(seed)
            Random.seed!(seed)
        end
        ds = Index.(ntuple(x -> d, 2))
        χs = Index.(ntuple(x -> χ, 4))
        A = randomITensor(χs[1], ds[1], χs[2])
        A = A + swapind(A, χs[1], χs[2])
        At = copy(A)
        At = swapinds(At, (χs[1], ds[1], χs[2]), (χs[3], ds[2], χs[4]))
        T = At * A * delta(χs[2], χs[3])
        Di = combiner(ds[1],ds[2]; tags="Di")
        Χi = combiner(χs[1],χs[4]; tags="Χi")
        di = combinedind(Di)
        χi = combinedind(Χi)
        E = T * Di * Χi
        U, S, V = svd(E, (di))

        Up = swapprime(dag(U), 0 => 1) * U * delta(inds(S)[1], inds(S)[1]')
        P = Up - delta(inds(Up))
        return P |> matrix
    end
end
