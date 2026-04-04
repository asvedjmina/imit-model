using DrWatson
@quickactivate "lab_04_models"

ENV["GKSwstype"] = "100"

include(srcdir("sir_analysis.jl"))

params = default_parameters()
agent_df, model_df = run_basic_experiment(params)
basic_plot = plot_basic_dynamics(agent_df, model_df)
basic_plot

savefig(basic_plot, plotsdir("sir_basic_dynamics.png"))
@save datadir("sir_basic_agent.jld2") agent_df
@save datadir("sir_basic_model.jld2") model_df

beta_df, beta_grouped = scan_beta_experiment()
beta_plot = plot_beta_scan(beta_grouped)
beta_plot

CSV.write(datadir("beta_scan_all.csv"), beta_df)
savefig(beta_plot, plotsdir("beta_scan.png"))

migration_df, migration_grouped = scan_migration_experiment()
migration_plot = plot_migration_effect(migration_grouped)
migration_plot

CSV.write(datadir("migration_scan_all.csv"), migration_df)
savefig(migration_plot, plotsdir("migration_effect.png"))

comprehensive_plot = comprehensive_analysis_plot(beta_df)
comprehensive_plot

savefig(comprehensive_plot, plotsdir("comprehensive_analysis.png"))
