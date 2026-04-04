using DrWatson
@quickactivate "lab_04_models"

ENV["GKSwstype"] = "100"

include(srcdir("sir_analysis.jl"))

base_global, base_city, hetero_global, hetero_city, summary_df, city_summary =
    heterogeneity_experiment()

summary_df

global_plot = plot_heterogeneity_global(base_global, hetero_global)
global_plot

savefig(global_plot, plotsdir("heterogeneity_global_comparison.png"))

city_plot = plot_city_dynamics(hetero_city; title = "Неоднородный beta по городам")
city_plot

savefig(city_plot, plotsdir("heterogeneity_city_dynamics.png"))

combined_plot = plot(global_plot, city_plot; layout = (2, 1), size = (900, 800))
combined_plot

CSV.write(datadir("heterogeneity_summary.csv"), summary_df)
CSV.write(datadir("heterogeneity_city_summary.csv"), city_summary)
savefig(combined_plot, plotsdir("heterogeneity_effect.png"))
