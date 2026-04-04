using DrWatson
@quickactivate "lab_04_models"

ENV["GKSwstype"] = "100"

include(srcdir("sir_analysis.jl"))

migration_df, migration_grouped = scan_migration_experiment()
migration_plot = plot_migration_effect(migration_grouped)

CSV.write(datadir("migration_scan_all.csv"), migration_df)
savefig(migration_plot, plotsdir("migration_effect.png"))

println("Migration scan completed.")
