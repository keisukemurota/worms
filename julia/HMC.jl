using LinearAlgebra
using Plots
using Zygote

const μ = 0.0
const σ = 1.0

function logprob(x)
	return -(x - 1)^2 * (x + 1)^2
end

plot(x -> logprob(x) |> exp, -10, 10)

function dlogprob(x)
	return Zygote.gradient(x -> logprob(x), x)[1]
	# return -(x-μ) / (σ^2)
end

function momentum(p, τ)
	return (p .^ 2) / (2 * τ^2)
end

function dmomentum(p, τ)
	return p ./ (τ^2)
end

function H(p::point, τ)
	return -logprob(p.x) + momentum(p.p, τ)
end

mutable struct point
	x::Float64
	p::Float64
end


function leapfrog!(p::point, ϵ, τ)
	p.p -= -0.5ϵ * dlogprob(p.x)
	p.x += ϵ * dmomentum(p.p, τ)
	p.p -= -0.5ϵ * dlogprob(p.x)
end

function HMC(x, τ, ϵ, T)::Float64
	p = point(x, randn())
	h = H(p, τ)
	for i in 1:T
		leapfrog!(p, ϵ, τ)
	end
	r = rand()
	if r < exp(-H(p, τ) + h)
		return p.x
	else
		return x
	end
end


ϵ = 0.01
τ = 1

# steps = 100
# anim = @animate for i in 1:steps
#     leapfrog!(p, ϵ, τ)
#     println(H(p, τ))
#     plot([p.x], [p.p], seriestype=:scatter, xlims=(-10, 10), ylims=(-10, 10), title="Position at Step $i")
# end
# gif(anim, "leapfrog_animation.gif", fps = 15)


x = 0.0
xs = Vector{Float64}([x])
for i in 2:100000
	x = HMC(x, 1, 0.1, 30)
	push!(xs, x)
end

using StatsPlots

density(xs)
plot!(x -> logprob(x) |> exp, color = :red)
