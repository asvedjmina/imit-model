using DataFrames
using DrWatson
using Random
using Test

@quickactivate "lab_06_models"

include(srcdir("SIRPetri.jl"))
using .SIRPetri

@testset "lab_06_models tests" begin
    net, u0, states = build_sir_network()
    @test states == [:S, :I, :R]
    @test length(u0) == 3
    @test sum(u0) == 1000.0

    df_det = simulate_deterministic(net, u0, (0.0, 20.0); saveat = 0.5)
    @test nrow(df_det) > 10
    @test df_det.time[1] == 0.0
    @test df_det.time[end] == 20.0
    @test all(df_det.S .>= 0.0)
    @test abs((df_det.S[end] + df_det.I[end] + df_det.R[end]) - 1000.0) < 1e-6

    df_stoch = simulate_stochastic(net, u0, (0.0, 20.0); rng = MersenneTwister(42))
    @test nrow(df_stoch) >= 2
    @test df_stoch.time[1] == 0.0
    @test df_stoch.time[end] == 20.0
    @test all(df_stoch.I .>= 0.0)

    summary = summarize_trajectory(df_det)
    @test summary.peak_I >= 10.0
    @test summary.final_R > 0.0

    scan_df = parameter_scan(0.10:0.10:0.30; gamma = 0.1, tmax = 40.0, saveat = 1.0)
    @test nrow(scan_df) == 3
    @test :peak_I in Symbol.(names(scan_df))

    grid_df = parameter_grid([(beta = 0.2, gamma = 0.1), (beta = 0.4, gamma = 0.1)])
    @test nrow(grid_df) == 2
    @test grid_df.peak_I[2] > grid_df.peak_I[1]

    dot = to_graphviz_sir(net)
    @test occursin("infection", dot)
    @test occursin("recovery", dot)
end
