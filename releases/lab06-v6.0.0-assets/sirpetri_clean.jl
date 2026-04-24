using DrWatson
@quickactivate "lab_06_models"

ENV["GKSwstype"] = "100"

include(joinpath(@__DIR__, "..", "src", "SIRPetri.jl"))
using .SIRPetri
using CSV
using DataFrames
using Plots
using Random

mkpath(datadir())
mkpath(plotsdir())

beta = 0.30
gamma = 0.10
tmax = 100.0
saveat = 0.5
seed = 123

net, u0, states = build_sir_network(beta, gamma)

network_plot = plot_sir_network(net)
network_plot

savefig(network_plot, plotsdir("sir_network.png"))
write(joinpath(plotsdir(), "sir_network.dot"), to_graphviz_sir(net))

df_det = simulate_deterministic(
    net,
    u0,
    (0.0, tmax);
    saveat = saveat,
    rates = [beta, gamma],
)

det_summary = DataFrame([summarize_trajectory(df_det)])
det_summary

CSV.write(datadir("sir_det.csv"), df_det)
CSV.write(datadir("sir_det_summary.csv"), det_summary)

det_plot = plot_sir(df_det; title = "Deterministic SIR dynamics")
det_plot

savefig(det_plot, plotsdir("sir_det_dynamics.png"))

Random.seed!(seed)
df_stoch = simulate_stochastic(
    net,
    u0,
    (0.0, tmax);
    rates = [beta, gamma],
    rng = MersenneTwister(seed),
)

stoch_summary = DataFrame([summarize_trajectory(df_stoch)])
stoch_summary

CSV.write(datadir("sir_stoch.csv"), df_stoch)
CSV.write(datadir("sir_stoch_summary.csv"), stoch_summary)

stoch_plot = plot_sir(df_stoch; title = "Stochastic SIR dynamics")
stoch_plot

savefig(stoch_plot, plotsdir("sir_stoch_dynamics.png"))

comparison_plot = plot_comparison(df_det, df_stoch)
comparison_plot

savefig(comparison_plot, plotsdir("comparison.png"))

beta_range = 0.10:0.05:0.80
df_scan = parameter_scan(beta_range; gamma = gamma, tmax = tmax, saveat = saveat)
df_scan

CSV.write(datadir("sir_scan.csv"), df_scan)

scan_plot = plot_scan(df_scan)
scan_plot

savefig(scan_plot, plotsdir("sir_scan.png"))

sensitivity_plot = plot_sensitivity(df_scan)
sensitivity_plot

savefig(sensitivity_plot, plotsdir("sensitivity.png"))

df_anim = simulate_deterministic(
    net,
    u0,
    (0.0, tmax);
    saveat = 1.0,
    rates = [beta, gamma],
)

animate_sir(df_anim, plotsdir("sir_animation.gif"); fps = 8)
