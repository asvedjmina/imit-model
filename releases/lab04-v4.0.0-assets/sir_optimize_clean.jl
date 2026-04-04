using DrWatson
@quickactivate "lab_04_models"

ENV["GKSwstype"] = "100"

include(srcdir("sir_analysis.jl"))

peak_limit = 0.30
all_df, feasible_df, best = optimize_under_peak_constraint(; peak_limit)

first(feasible_df, min(10, nrow(feasible_df)))

opt_plot = plot_constrained_optimization(all_df, feasible_df; peak_limit)
opt_plot

@save datadir("optimization_result.jld2") best feasible_df all_df
CSV.write(datadir("optimization_samples.csv"), all_df)
CSV.write(datadir("optimization_pareto.csv"), feasible_df)
savefig(opt_plot, plotsdir("optimization_pareto.png"))
