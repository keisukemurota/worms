using JLD2
using Plots

file = joinpath(dirname(@__FILE__),"pickles", "res.jld2")
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


(res_list |> keys |> collect)[148]

(res_list |> values |> collect)[148].spec_sys.losses_track

η1


scatter(η1, η2, 
    xlabel="orthogonal", 
    ylabel="special linear", 
    title="Negativity optimized with orthogonal and special linear",
    size=(1000, 1000),  # Increase plot size to 800x600 pixels
    markersize=5,  # Optionally increase marker size for better visibility
    fontsize=40,  # Increase font size for better readability
    titlefontsize=20,  # Increase title font size
    guidefontsize=20,  # Increase font size for axis labels
    tickfontsize=10,  # Increase font size for tick labels
    legend=false,
    aspect_ratio = 1,
    xlims = (0, 1.5),
    ylims = (0, 1.5),
)


scatter(η′1, η′2, 
    xlabel="orthogonal", 
    ylabel="unitary", 
    title="Adaptive L1 loss optimized with orthogonal and unitary",
    size=(800, 600),  # Increase plot size to 800x600 pixels
    markersize=5,  # Optionally increase marker size for better visibility
    fontsize=40,  # Increase font size for better readability
    titlefontsize=16,  # Increase title font size
    guidefontsize=16,  # Increase font size for axis labels
    tickfontsize=10,  # Increase font size for tick labels
    legend=false,
)

file = joinpath(dirname(@__FILE__),"pickles", "spec_vs_uni.jld2")
# current path

res_list = load_object(file);

z = Float64[]
function color(jx, jy)
    r = res[(jx, jy)]
    if r.spec_sys == []
        return 0.0
    else 
        # tmp = - (r.spec_sys.losses_opt |> minimum) + (r.orth_sys.losses_opt |> minimum)
        # tmp = r.spec_sys.losses_opt |> minimum
        tmp = r.uni_sys.losses_opt |> minimum
    end
    return tmp
end

color(0.0, 0.0)

Jy = 0:0.1:2
Jz = 0:0.1:2
# using Plots
heatmap(
    Jy, 
    Jz, 
    [color(jz, jy) for jz in Jz, jy in Jy], 
    size = (700, 700),
    aspect_ratio = :equal,
    xlims = (0, 2),
    ylims = (0, 2),
    xlabel = "Jz",
    ylabel = "Jy",
    title = "L1 adaptive loss (unitary)",
    colorbar = true,
    c = :Reds,
    # margin = 0Plots.mm,
    rightmargin = 2Plots.mm,
    titlefontsize = 20,
    guidefontsize = 12,
    tickfontsize = 12,
    legendfontsize = 12,
    clims = (0, 1.0),
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