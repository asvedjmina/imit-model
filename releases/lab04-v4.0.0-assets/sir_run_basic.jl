using DrWatson
@quickactivate "lab_04_models"

ENV["GKSwstype"] = "100"

include(srcdir("sir_analysis.jl"))

params = default_parameters()
agent_df, model_df = run_basic_experiment(params)
basic_plot = plot_basic_dynamics(agent_df, model_df)
summary_df = DataFrame([basic_summary(params)])

savefig(basic_plot, plotsdir("sir_basic_dynamics.png"))
@save datadir("sir_basic_agent.jld2") agent_df
@save datadir("sir_basic_model.jld2") model_df
CSV.write(datadir("basic_summary.csv"), summary_df)

println("Basic experiment completed.")
