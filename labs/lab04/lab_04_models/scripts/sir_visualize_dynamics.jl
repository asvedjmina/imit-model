using DrWatson
@quickactivate "lab_04_models"

ENV["GKSwstype"] = "100"

include(srcdir("sir_analysis.jl"))

beta_df = CSV.read(datadir("beta_scan_all.csv"), DataFrame)
comprehensive_plot = comprehensive_analysis_plot(beta_df)

savefig(comprehensive_plot, plotsdir("comprehensive_analysis.png"))

println("Comprehensive plot completed.")
