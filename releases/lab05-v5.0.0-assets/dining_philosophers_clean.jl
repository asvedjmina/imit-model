using DrWatson
@quickactivate "lab_05_models"

ENV["GKSwstype"] = "100"

include(joinpath(@__DIR__, "..", "src", "DiningPhilosophers.jl"))
using .DiningPhilosophers
using CSV
using DataFrames
using Plots
using Random

mkpath(datadir())
mkpath(plotsdir())

N = 5
tmax = 50.0
classic_seed = 2026
arbiter_seed = 2027

net_classic, u0_classic, _ = build_classical_network(N)
net_arbiter, u0_arbiter, _ = build_arbiter_network(N)

classic_network_plot = plot_network(net_classic)
arbiter_network_plot = plot_network(net_arbiter)

classic_network_plot
arbiter_network_plot

savefig(classic_network_plot, plotsdir("classic_network.png"))
savefig(arbiter_network_plot, plotsdir("arbiter_network.png"))

df_classic = simulate_stochastic(
    net_classic,
    u0_classic,
    tmax;
    rng = MersenneTwister(classic_seed),
)
df_arbiter = simulate_stochastic(
    net_arbiter,
    u0_arbiter,
    tmax;
    rng = MersenneTwister(arbiter_seed),
)

CSV.write(datadir("dining_classic.csv"), df_classic)
CSV.write(datadir("dining_arbiter.csv"), df_arbiter)

deadlock_classic = detect_deadlock(df_classic, net_classic)
deadlock_arbiter = detect_deadlock(df_arbiter, net_arbiter)

deadlock_summary = DataFrame(
    model = ["classic", "arbiter"],
    deadlock = [deadlock_classic, deadlock_arbiter],
    deadlock_time = [
        deadlock_time(df_classic, net_classic),
        deadlock_time(df_arbiter, net_arbiter),
    ],
    final_total_eaters = [
        total_eaters(df_classic, N)[end],
        total_eaters(df_arbiter, N)[end],
    ],
)
deadlock_summary

CSV.write(datadir("deadlock_summary.csv"), deadlock_summary)

classic_marking_plot = plot_marking_evolution(df_classic, net_classic)
arbiter_marking_plot = plot_marking_evolution(df_arbiter, net_arbiter)

classic_marking_plot
arbiter_marking_plot

savefig(classic_marking_plot, plotsdir("classic_simulation.png"))
savefig(arbiter_marking_plot, plotsdir("arbiter_simulation.png"))

df_classic_det = simulate_deterministic(net_classic, u0_classic, tmax)
df_arbiter_det = simulate_deterministic(net_arbiter, u0_arbiter, tmax)

deterministic_plot = plot_total_eaters_comparison(
    df_classic,
    df_classic_det,
    df_arbiter,
    df_arbiter_det,
    N,
)
deterministic_plot

CSV.write(datadir("dining_classic_deterministic.csv"), df_classic_det)
CSV.write(datadir("dining_arbiter_deterministic.csv"), df_arbiter_det)
savefig(deterministic_plot, plotsdir("deterministic_comparison.png"))

final_report_plot = plot_eating_comparison(df_classic, df_arbiter, N)
final_report_plot

savefig(final_report_plot, plotsdir("final_report.png"))

net_anim, u0_anim, _ = build_classical_network(3)
df_anim = simulate_stochastic(net_anim, u0_anim, 30.0; rng = MersenneTwister(123))

animate_marking(
    df_anim,
    net_anim,
    plotsdir("philosophers_simulation.gif");
    fps = 2,
)
