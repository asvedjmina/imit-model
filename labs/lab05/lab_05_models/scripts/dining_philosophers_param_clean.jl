using DrWatson
@quickactivate "lab_05_models"

ENV["GKSwstype"] = "100"

include(joinpath(@__DIR__, "..", "src", "DiningPhilosophers.jl"))
using .DiningPhilosophers
using CSV
using DataFrames
using Plots

mkpath(datadir())
mkpath(plotsdir())

philosopher_counts = 3:7
seeds = 1:16
tmax = 60.0

raw_scan = parameter_scan(philosopher_counts; seeds = seeds, tmax = tmax)
summary_scan = summarize_parameter_scan(raw_scan)

raw_scan
summary_scan

CSV.write(datadir("parameter_scan_raw.csv"), raw_scan)
CSV.write(datadir("parameter_scan_summary.csv"), summary_scan)

parameter_plot = plot_parameter_summary(summary_scan)
parameter_plot

savefig(parameter_plot, plotsdir("parameter_summary.png"))
