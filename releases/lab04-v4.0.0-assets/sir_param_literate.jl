# # Parameter scenarios for the graph SIR model
#
# This script extends the base experiment with several parameter scenarios:
# heterogeneous transmission, faster migration, and slower detection.
#
# ## Libraries and project activation

using DrWatson
@quickactivate "lab_04_models"

ENV["GKSwstype"] = "100"

include(srcdir("sir_analysis.jl"))

# ## Scenario list
#
# Each scenario modifies a small set of parameters while keeping the model
# structure unchanged.

scenarios = scenario_parameters()
summaries = NamedTuple[]

for spec in scenarios
    global_df, city_df, summary = run_scenario(spec)
    scenario_plot = plot_scenario(global_df, city_df; title = spec[:label])
    scenario_plot
    plot_name = "scenario_" * spec[:label] * ".png"

#jl     savefig(scenario_plot, plotsdir(plot_name))

    push!(summaries, summary)
end

summary_df = DataFrame(summaries)
summary_df

#jl CSV.write(datadir("scenario_summary.csv"), summary_df)

# ## Scenario comparison
#
# We compare peak infection fractions and death counts across all scenarios.

comparison_plot = scenario_comparison_plot(summary_df)
comparison_plot

#jl savefig(comparison_plot, plotsdir("scenario_comparison.png"))

# ## Conclusions
#
# 1. Heterogeneous transmission shifts the infection burden toward the first city.
# 2. Faster migration reduces the time to epidemic peak.
# 3. Slower detection increases both the peak load and the number of deaths.
