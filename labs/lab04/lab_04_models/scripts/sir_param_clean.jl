using DrWatson
@quickactivate "lab_04_models"

ENV["GKSwstype"] = "100"

include(srcdir("sir_analysis.jl"))

scenarios = scenario_parameters()
summaries = NamedTuple[]

for spec in scenarios
    global_df, city_df, summary = run_scenario(spec)
    scenario_plot = plot_scenario(global_df, city_df; title = spec[:label])
    scenario_plot
    plot_name = "scenario_" * spec[:label] * ".png"

    savefig(scenario_plot, plotsdir(plot_name))

    push!(summaries, summary)
end

summary_df = DataFrame(summaries)
summary_df

CSV.write(datadir("scenario_summary.csv"), summary_df)

comparison_plot = scenario_comparison_plot(summary_df)
comparison_plot

savefig(comparison_plot, plotsdir("scenario_comparison.png"))
