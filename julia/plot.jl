using Plots
using LinearAlgebra

function extract_data!(means::Vector{Float64}, errors::Vector{Float64}, line::String, keyword::String)
	if contains(line, keyword)
		# Extract the part after the keyword and before the next comma
		data_part = split(split(line, keyword)[2], ",")[1]
		# Remove any leading/trailing whitespace and the equals sign
		data_part = strip(data_part, [' ', '='])
		# Split the string on " +- " to separate mean and error
		parts = split(data_part, " +- ")
		if length(parts) == 2
			mean_val, error_val = parts
		else
			mean_val = parts[1]
			error_val = "0.0"
		end
		# Parse the cleaned strings as Float64 and push to the respective vectors
		push!(means, parse(Float64, mean_val))
		push!(errors, parse(Float64, error_val))
	end
end

function read_data(filename::String)
	Δc_list = Float64[]
	c_list = Float64[]
	Δe_list = Float64[]
	e_list = Float64[]
	as_list = Float64[]
	Δas_list = Float64[]
	T_list = Float64[]
	ΔT_list = Float64[]

	open(filename) do file
		for line in eachline(file)
			extract_data!(c_list, Δc_list, line, "Specific heat")
			extract_data!(e_list, Δe_list, line, "Energy per site")
			extract_data!(as_list, Δas_list, line, "Average sign")
			extract_data!(T_list, ΔT_list, line, "T")
		end
	end

	return T_list, c_list, Δc_list, e_list, Δe_list, as_list, Δas_list
end

T_list_optim, c_list_optim, Δc_list_optim, e_list_optim, Δe_list_optim, as_list_optim, Δas_list_optim = read_data("python/visualize/optim.txt")
T_list_none, c_list_none, Δc_list_none, e_list_none, Δe_list_none, as_list_none, Δas_list_none = read_data("python/visualize/none.txt")

p1 = plot(
	T_list_optim,
	c_list_optim,
	yerror = Δc_list_optim,
	marker = :circle,
	markersize = 3,
	ylim = (-0.1, 0.6),
	xlabel = "Temperature (T)",
	ylabel = "Specific Heat",
	label = "optim",
    markerstrokecolor=:auto
)
plot!(
	T_list_none,
	c_list_none,
	yerror = Δc_list_none,
	marker = :circle,
	markersize = 3,
	label = "none",
    markerstrokecolor=:auto
)

p2 = plot(
	T_list_optim,
	e_list_optim,
	yerror = Δe_list_optim,
	marker = :circle,
	markersize = 3,
	ylim = (-0.2, 1.2),
	xlabel = "Temperature (T)",
	ylabel = "Energy per Site",
	label = "optim",
    markerstrokecolor=:auto
)
plot!(
	T_list_none,
	e_list_none,
	yerror = Δe_list_none,
	marker = :circle,
	markersize = 3,
	label = "none",
    markerstrokecolor=:auto
)

p3 = plot(
	T_list_optim,
	as_list_optim,
	yerror = Δas_list_optim,
	marker = :circle,
	markersize = 3,
	xlabel = "Temperature (T)",
	ylabel = "Average Sign",
	ylim = (0, 1.2),
	label = "optim",
markerstrokecolor=:auto
)
plot!(
	T_list_none,
	as_list_none,
	yerror = Δas_list_none,
	marker = :circle,
	markersize = 3,
	label = "none",
    markerstrokecolor=:auto
)

display(p1)
display(p2)
display(p3)

