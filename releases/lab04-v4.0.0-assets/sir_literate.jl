# # SIR model on a graph of cities
#
# This literate script reproduces the main experiments for laboratory work 4:
# a base SIR run, a scan over transmission rate, a migration scan,
# and a summary visualization.
#
# ## Libraries and project activation

using DrWatson
@quickactivate "lab_04_models"

ENV["GKSwstype"] = "100"

include(srcdir("sir_analysis.jl"))

# ## Base experiment
#
# We start with three cities of equal population and seed one infected agent
# in the third city.

params = default_parameters()
agent_df, model_df = run_basic_experiment(params)
basic_plot = plot_basic_dynamics(agent_df, model_df)
basic_plot

#jl savefig(basic_plot, plotsdir("sir_basic_dynamics.png"))
#jl @save datadir("sir_basic_agent.jld2") agent_df
#jl @save datadir("sir_basic_model.jld2") model_df

# ## Scan over transmission intensity
#
# The next experiment varies `beta` from `0.1` to `1.0`
# and averages the result over three random seeds.

beta_df, beta_grouped = scan_beta_experiment()
beta_plot = plot_beta_scan(beta_grouped)
beta_plot

#jl CSV.write(datadir("beta_scan_all.csv"), beta_df)
#jl savefig(beta_plot, plotsdir("beta_scan.png"))

# ## Migration effect
#
# We vary the migration intensity and observe how it changes the time to peak
# and the peak number of infected agents.

migration_df, migration_grouped = scan_migration_experiment()
migration_plot = plot_migration_effect(migration_grouped)
migration_plot

#jl CSV.write(datadir("migration_scan_all.csv"), migration_df)
#jl savefig(migration_plot, plotsdir("migration_effect.png"))

# ## Summary plot
#
# The final figure combines three panels built from the beta scan:
# peak infection, deaths, and the final recovered fraction.

comprehensive_plot = comprehensive_analysis_plot(beta_df)
comprehensive_plot

#jl savefig(comprehensive_plot, plotsdir("comprehensive_analysis.png"))

# ## Conclusions
#
# 1. The base run shows a standard epidemic wave with a clear infection peak.
# 2. The threshold behavior becomes visible in the beta scan.
# 3. Higher migration moves the infection peak earlier and makes it stronger.
# 4. The combined figure is suitable for the lab report and presentation.
