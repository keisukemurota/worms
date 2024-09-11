using LinearAlgebra


transition_matrix = randn(5, 5)

transition_matrix = abs.(transition_matrix)


transition_matrix ./= sum(transition_matrix, dims=1)

transition_matrix = transition_matrix'

transition_matrix[1, :] |> sum

tm = transition_matrix

cumsum(tm, dims=2)

tm[3, :] .= 0
tm[3, 4] = 1

tm[3, :] |> sum

v = tm[1, :]

i = 

