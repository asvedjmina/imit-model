using DrWatson
@quickactivate "lab_04_models"

ENV["GKSwstype"] = "100"

include(srcdir("sir_analysis.jl"))

peak_limit = 0.30
all_df, feasible_df, best = optimize_under_peak_constraint(; peak_limit)
opt_plot = plot_constrained_optimization(all_df, feasible_df; peak_limit)

@save datadir("optimization_result.jld2") best feasible_df all_df
CSV.write(datadir("optimization_samples.csv"), all_df)
CSV.write(datadir("optimization_pareto.csv"), feasible_df)
savefig(opt_plot, plotsdir("optimization_pareto.png"))

println("Optimization completed.")
if best === nothing
    println("No feasible solution found under peak limit $(peak_limit).")
else
    println("Best feasible candidate:")
    println("  beta_und            = $(round(best.beta_und, digits = 4))")
    println("  detection_time      = $(best.detection_time)")
    println("  death_rate          = $(round(best.death_rate, digits = 4))")
    println("  mean_peak           = $(round(best.mean_peak, digits = 4))")
    println("  mean_deaths         = $(round(best.mean_deaths, digits = 2))")
end
