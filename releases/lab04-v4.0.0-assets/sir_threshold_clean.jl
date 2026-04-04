using DrWatson
@quickactivate "lab_04_models"

ENV["GKSwstype"] = "100"

include(srcdir("sir_analysis.jl"))

threshold_df, threshold_grouped, threshold_summary = threshold_study()
threshold_grouped

threshold_plot = plot_threshold_study(
    threshold_grouped;
    epidemic_threshold = threshold_summary.epidemic_threshold,
    theoretical_beta = threshold_summary.theoretical_beta,
    observed_beta = threshold_summary.observed_beta,
)
threshold_plot

CSV.write(datadir("threshold_scan_all.csv"), threshold_df)
CSV.write(datadir("threshold_scan_grouped.csv"), threshold_grouped)
CSV.write(datadir("threshold_summary.csv"), DataFrame([threshold_summary]))
savefig(threshold_plot, plotsdir("threshold_study.png"))
