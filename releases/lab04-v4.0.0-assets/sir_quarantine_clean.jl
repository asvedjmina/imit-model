using DrWatson
@quickactivate "lab_04_models"

ENV["GKSwstype"] = "100"

include(srcdir("sir_analysis.jl"))

summary_df, dynamics = quarantine_experiment()
summary_df

summary_plot = plot_quarantine_summary(summary_df)
summary_plot

CSV.write(datadir("quarantine_summary.csv"), summary_df)
savefig(summary_plot, plotsdir("quarantine_effect.png"))

candidate_df = filter(row -> row.scenario != "no_quarantine" && row.locked_cities > 0, summary_df)
best_label = if nrow(candidate_df) == 0
    "no_quarantine"
else
    candidate_df[argmin(candidate_df.deaths), :scenario]
end

baseline_df = dynamics["no_quarantine"].global_df
best_df = dynamics[best_label].global_df
dynamics_plot = plot_quarantine_dynamics(baseline_df, best_df; quarantine_label = best_label)
dynamics_plot

savefig(dynamics_plot, plotsdir("quarantine_best_dynamics.png"))
