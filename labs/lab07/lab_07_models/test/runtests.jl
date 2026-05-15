using DataFrames
using DrWatson
using Random
using Test

@quickactivate "lab_07_models"

include(srcdir("DiscreteEventModels.jl"))
using .DiscreteEventModels

@testset "lab_07_models tests" begin
    mmc_params = MMcParameters(lambda = 0.45, mu = 0.5, servers = 2, num_customers = 250)
    analytic = erlang_c_metrics(mmc_params)

    @test isapprox(analytic.rho, 0.45; atol = 1.0e-12)
    @test 0.0 < analytic.p_wait < 1.0
    @test analytic.W > analytic.Wq > 0.0

    mmc_df = simulate_mmc(mmc_params; rng = MersenneTwister(42))
    mmc_summary = summarize_mmc(mmc_df, mmc_params)
    mmc_comparison = compare_mmc_metrics(mmc_df, mmc_params)

    @test nrow(mmc_df) == mmc_params.num_customers
    @test all(mmc_df.waiting_time .>= 0.0)
    @test all(mmc_df.departure_time .>= mmc_df.service_start .>= mmc_df.arrival_time)
    @test maximum(mmc_df.server_id) <= mmc_params.servers
    @test 0.0 <= mmc_summary.wait_probability <= 1.0
    @test nrow(mmc_comparison) == 5
    @test abs(mmc_summary.avg_waiting_time - analytic.Wq) < 0.8

    mmc_sweep = mmc_server_sweep(1:3; lambda = 0.45, mu = 0.5, num_customers = 200, seed = 11)
    @test nrow(mmc_sweep) == 3
    @test mmc_sweep.analytic_wait_probability[1] > mmc_sweep.analytic_wait_probability[end]

    ross_base = RossParameters(machines = 10, spares = 3, repairers = 1, failure_mean = 100.0, repair_mean = 1.0)
    ross_more_spares = RossParameters(machines = 10, spares = 4, repairers = 1, failure_mean = 100.0, repair_mean = 1.0)
    ross_more_repairers = RossParameters(machines = 10, spares = 3, repairers = 2, failure_mean = 100.0, repair_mean = 1.0)

    analytic_base = ross_mean_time_to_failure_analytic(ross_base)
    analytic_more_spares = ross_mean_time_to_failure_analytic(ross_more_spares)
    analytic_more_repairers = ross_mean_time_to_failure_analytic(ross_more_repairers)

    @test analytic_base > 0.0
    @test analytic_more_spares > analytic_base
    @test analytic_more_repairers > analytic_base

    ross_result = simulate_ross(ross_base; seed = 17)
    @test ross_result.crash_time > 0.0
    @test 0.0 <= ross_result.utilization <= 1.0
    @test ross_result.mean_queue_length >= 0.0
    @test ross_result.trace.healthy_machines[1] == ross_base.machines + ross_base.spares
    @test ross_result.trace.healthy_machines[end] == ross_base.machines - 1
    @test ross_result.trace.event[end] == "crash"

    ross_runs = ross_replicates(ross_base; runs = 4, seed = 200)
    @test nrow(ross_runs) == 4
    @test all(ross_runs.crash_time .> 0.0)

    machine_sweep = ross_machine_sweep([8, 10]; spares = 3, repairers = 2, runs = 4, seed = 300)
    @test nrow(machine_sweep) == 2
    @test machine_sweep.analytic_crash_time[1] > machine_sweep.analytic_crash_time[2]

    repairer_sweep = ross_repairer_sweep(1:2; machines = 10, spares = 3, runs = 4, seed = 400)
    @test nrow(repairer_sweep) == 2
    @test repairer_sweep.analytic_crash_time[2] > repairer_sweep.analytic_crash_time[1]
end
