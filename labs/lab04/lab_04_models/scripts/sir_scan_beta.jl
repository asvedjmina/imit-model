using DrWatson
@quickactivate "lab_04_models"

ENV["GKSwstype"] = "100"

include(srcdir("sir_analysis.jl"))

beta_df, beta_grouped = scan_beta_experiment()
beta_plot = plot_beta_scan(beta_grouped)

CSV.write(datadir("beta_scan_all.csv"), beta_df)
savefig(beta_plot, plotsdir("beta_scan.png"))

println("Beta scan completed.")
