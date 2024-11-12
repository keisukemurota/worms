using JLD2
using Plots
using LinearAlgebra
Pkg.activate("./MinEigLoss/MinEigLoss")
using MinEigLoss

file = joinpath(dirname(@__FILE__), "pickles", "res.jld2")
# current path

res_list = load_object(file);

begin
	η1 = []
	η2 = []
	η′1 = []
	η′2 = []
	for (key, res) in res_list
		tmp1 = res.orth_sys.losses_opt |> minimum
		tmp2 = res.spec_sys.losses_opt |> minimum
		if tmp2 < tmp1
			println("key: $key, tmp1: $tmp1, tmp2: $tmp2")
		end
		push!(η1, tmp1)
		push!(η2, min(tmp1, tmp2))
		push!(η′1, res.orth_loc.losses_opt |> minimum)
		push!(η′2, res.uni_loc.losses_opt |> minimum)
	end
end

η1 - η2 |> maximum


(res_list|>keys|>collect)[148]

(res_list|>values|>collect)[148].spec_sys.losses_track

η1

function stochastic(A)
    B = copy(A)
    for i in 1:size(A, 1)
        for j in 1:size(A, 2)
            if i != j
                B[i, j] = -abs(A[i, j])
            end
        end
    end
    return B
end

(U * res_list[-1.0, -1.0].H * U')  |> stochastic |> eigmax
res_list[-1.0, -1.0].H  |> eigmax
H = res_list[-1.0, -1.0].H

u = res_list[-1.0, -1.0].orth_sys.thetas[1]
U = kron(u, u, u)

res_list[-1.0, -1.0].orth_sys.losses_opt |> minimum

loss = MinEigLoss.mle_sys_unitary(H - 100I(8), 2)
loss(I(2))

color(-1.0, -1.0, type = :orth)


scatter(η1, η2,
	xlabel = "orthogonal",
	ylabel = "special linear",
	title = "Negativity optimized with orthogonal and special linear",
	size = (1000, 1000),  # Increase plot size to 800x600 pixels
	markersize = 5,  # Optionally increase marker size for better visibility
	fontsize = 40,  # Increase font size for better readability
	titlefontsize = 20,  # Increase title font size
	guidefontsize = 20,  # Increase font size for axis labels
	tickfontsize = 10,  # Increase font size for tick labels
	legend = false,
	aspect_ratio = 1,
	xlims = (0, 1.5),
	ylims = (0, 1.5),
)


scatter(η′1, η′2,
	xlabel = "orthogonal",
	ylabel = "unitary",
	title = "Adaptive L1 loss optimized with orthogonal and unitary",
	size = (800, 600),  # Increase plot size to 800x600 pixels
	markersize = 5,  # Optionally increase marker size for better visibility
	fontsize = 40,  # Increase font size for better readability
	titlefontsize = 16,  # Increase title font size
	guidefontsize = 16,  # Increase font size for axis labels
	tickfontsize = 10,  # Increase font size for tick labels
	legend = false,
)

file = joinpath(dirname(@__FILE__), "pickles", "spec_vs_uni_complex.jld2")
# current path

# res_list = load_object(file);
res_list = load_object("/Users/keisuke/Documents/projects/todo/worms/spec_vs_uni_complex.jld2");


res_list |>keys |> collect |> sort

res_list[2.0, 2.0].spec_sys_uni.losses_opt |> minimum
res_list[2.0, 2.0].spec_sys_orth.losses_opt |> minimum
res_list[2.0, 2.0].orth_sys.losses_opt |> minimum
res_list[2.0, 2.0].uni_sys.losses_opt |> minimum

z = Float64[]
function color(jx, jy; type = :orth)
	r = res_list[(jx, jy)]
	if r.orth_sys == []
		return 0.0
	else
		if type == :orth
			tmp = r.orth_sys.losses_opt |> minimum
		elseif type == :uni
			tmp = r.uni_sys.losses_opt |> minimum
		elseif type == :spec_orth
			tmp = r.spec_sys_orth.losses_opt |> minimum
		elseif type == :spec_uni
			tmp = r.spec_sys_uni.losses_opt |> minimum
		end
	end
	return tmp
end

color(0.0, 0.0, type = :spec_uni)

Jy = -2:0.1:2
Jz = -2:0.1:2
# using Plots

title_dict = Dict(:orth => "orthogonal", :uni => "unitary", :spec_orth => "special linear real", :spec_uni => "special linear complex")
begin
	plots = []
	for type in [:orth, :uni, :spec_orth, :spec_uni]
		p = heatmap(
			Jy,
			Jz,
			[color(jz, jy, type = type) for jz in Jz, jy in Jy],
			size = (600, 600),  # Reduced size
			aspect_ratio = :equal,
			xlims = (0, 2),
			ylims = (0, 2),
			xlabel = "Jx",
			ylabel = "Jy",
			title = title_dict[type],
			colorbar = false,
			# c = :Reds,
			# rightmargin = -4Plots.mm,  # Reduced right margin
			# leftmargin = 2Plots.mm,   # Increased left margin for y-label
			# bottommargin = -10Plots.mm, # Added bottom margin
			# topmargin = -20Plots.mm,    # Added top margin
			titlefontsize = 25,       # Reduced title font size
			guidefontsize = 15,       # Reduced guide font size
			tickfontsize = 13,         # Reduced tick font size
			legendfontsize = 10,       # Reduced legend font size
			clims = (0, 2.0),
			colorbar_titlefontsize = 10, # Set colorbar title font size
			framestyle = :box,        # Added frame around the plot
			yguide_position = :left,  # Ensure y-label is on the left
		)
		push!(plots, p)
	end
	plot(plots..., layout = (2, 2), size = (1400, 1400))
end


p = heatmap(
	Jy,
	Jz,
	[color(jz, jy, type = :uni) - color(jz, jy, type = :spec_uni) for jz in Jz, jy in Jy],
	size = (700, 700),
	aspect_ratio = :equal,
	xlims = (0, 2),
	ylims = (0, 2),
	xlabel = "Jz",
	ylabel = "Jy",
	title = "L1 adaptive loss (unitary - special linear complex)",
	colorbar = true,
	# c = :Reds,
	rightmargin = 2Plots.mm,
	titlefontsize = 20,
	guidefontsize = 12,
	tickfontsize = 12,
	legendfontsize = 12,
	# clims = (0, 1.0),
)

begin
	η1 = []
	η2 = []
	η′1 = []
	η′2 = []
	for (key, res) in res_list
		if res.orth_sys != []
			tmp1 = res.orth_sys.losses_opt |> minimum
			tmp2 = res.spec_sys.losses_opt |> minimum
			if tmp2 < tmp1
				println("key: $key, tmp1: $tmp1, tmp2: $tmp2")
			end
		else
			tmp1 = 0
			tmp2 = 0
		end
		push!(η1, tmp1)
		push!(η2, min(tmp1, tmp2))
	end
end

η1




res_list[0.0, 0.0].orth_sys