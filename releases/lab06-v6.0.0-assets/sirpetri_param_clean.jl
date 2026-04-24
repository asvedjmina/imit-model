using DrWatson
@quickactivate "lab_06_models"

ENV["GKSwstype"] = "100"

include(joinpath(@__DIR__, "..", "src", "SIRPetri.jl"))
using .SIRPetri
using CSV
using DataFrames
using Plots

mkpath(datadir())
mkpath(plotsdir())

parameter_pairs = [
    (beta = 0.15, gamma = 0.15),
    (beta = 0.20, gamma = 0.10),
    (beta = 0.25, gamma = 0.10),
    (beta = 0.30, gamma = 0.10),
    (beta = 0.35, gamma = 0.10),
    (beta = 0.40, gamma = 0.08),
    (beta = 0.45, gamma = 0.08),
    (beta = 0.50, gamma = 0.12),
    (beta = 0.55, gamma = 0.12),
]

grid_df = parameter_grid(parameter_pairs; tmax = 100.0, saveat = 0.5)
grid_df

CSV.write(datadir("sir_parameter_grid.csv"), grid_df)

peak_heatmap = plot_parameter_heatmap(grid_df; metric = :peak_I)
peak_heatmap

savefig(peak_heatmap, plotsdir("sir_parameter_heatmap.png"))

recovered_heatmap = plot_parameter_heatmap(grid_df; metric = :final_R)
recovered_heatmap

savefig(recovered_heatmap, plotsdir("sir_finalR_heatmap.png"))

sorted_df = sort(grid_df, :peak_I, rev = true)
labels = ["b=$(row.beta), g=$(row.gamma)" for row in eachrow(sorted_df)]

rank_plot = bar(
    labels,
    sorted_df.peak_I;
    xlabel = "Parameter pair",
    ylabel = "Peak I",
    title = "Peak infected for parameter pairs",
    legend = false,
    xrotation = 35,
    size = (900, 500),
)
rank_plot

savefig(rank_plot, plotsdir("sir_parameter_rank.png"))
