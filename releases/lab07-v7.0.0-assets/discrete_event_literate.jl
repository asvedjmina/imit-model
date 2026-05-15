# # Дискретно-событийные модели: M/M/c и модель Росса
#
# В лабораторной работе реализуются две системы:
#
# - многоканальная очередь `M/M/c`;
# - модель Росса с резервом и ремонтом машин.
#
# Код оформлен в структуре DrWatson и использует единый модуль
# для симуляции, аналитических расчётов и визуализации.

using DrWatson
@quickactivate "lab_07_models"

ENV["GKSwstype"] = "100"

include(joinpath(@__DIR__, "..", "src", "DiscreteEventModels.jl"))
using .DiscreteEventModels
using CSV
using DataFrames
using Plots
using Random
using Statistics

mkpath(datadir())
mkpath(plotsdir())

# ## Базовый эксперимент для M/M/c
#
# Используется система с двумя каналами обслуживания.
# Сначала выполняется событийная симуляция, затем результаты
# сравниваются с формулами Эрланга C.

mmc_params = MMcParameters(lambda = 0.9, mu = 0.5, servers = 2, num_customers = 800)
mmc_df = simulate_mmc(mmc_params; rng = MersenneTwister(123))
mmc_metrics = compare_mmc_metrics(mmc_df, mmc_params)
mmc_summary = DataFrame([merge(summarize_mmc(mmc_df, mmc_params), erlang_c_metrics(mmc_params))])

mmc_metrics
mmc_summary

#jl CSV.write(datadir("mmc_customers.csv"), mmc_df)
#jl CSV.write(datadir("mmc_metrics.csv"), mmc_metrics)
#jl CSV.write(datadir("mmc_summary.csv"), mmc_summary)

mmc_timeline = plot_mmc_timeline(mmc_df; title = "M/M/2: queue and service timeline")
mmc_timeline

#jl savefig(mmc_timeline, plotsdir("mmc_timeline.png"))

mmc_wait_hist = plot_mmc_wait_histogram(mmc_df; title = "M/M/2: waiting time histogram")
mmc_wait_hist

#jl savefig(mmc_wait_hist, plotsdir("mmc_wait_histogram.png"))

mmc_util_plot = plot_mmc_server_utilization(mmc_df, mmc_params; title = "M/M/2: server utilization")
mmc_util_plot

#jl savefig(mmc_util_plot, plotsdir("mmc_server_utilization.png"))

mmc_comparison = plot_mmc_comparison(mmc_df, mmc_params)
mmc_comparison

#jl savefig(mmc_comparison, plotsdir("mmc_sim_vs_analytic.png"))

# ## Базовый эксперимент для модели Росса
#
# Используется конфигурация с несколькими ремонтниками,
# чтобы выполнить расширение по условию лабораторной работы.

ross_params = RossParameters(machines = 10, spares = 3, repairers = 2, failure_mean = 100.0, repair_mean = 1.0)
ross_result = simulate_ross(ross_params; seed = 150)
ross_runs = ross_replicates(ross_params; runs = 40, seed = 500)
ross_repairers = ross_repairer_sweep(1:3; machines = 10, spares = 3, runs = 30, seed = 700)

ross_summary = DataFrame(
    [(
        machines = ross_params.machines,
        spares = ross_params.spares,
        repairers = ross_params.repairers,
        analytic_crash_time = ross_mean_time_to_failure_analytic(ross_params),
        sim_mean_crash_time = mean(ross_runs.crash_time),
        sim_std_crash_time = std(ross_runs.crash_time),
        mean_utilization = mean(ross_runs.utilization),
        mean_queue_length = mean(ross_runs.mean_queue_length),
    )],
)

ross_summary
ross_repairers

#jl CSV.write(datadir("ross_trace.csv"), ross_result.trace)
#jl CSV.write(datadir("ross_runs.csv"), ross_runs)
#jl CSV.write(datadir("ross_summary.csv"), ross_summary)
#jl CSV.write(datadir("ross_repairer_sweep.csv"), ross_repairers)

ross_trace_plot = plot_healthy_machines(ross_result.trace; title = "Ross model: healthy machines")
ross_trace_plot

#jl savefig(ross_trace_plot, plotsdir("ross_healthy_machines.png"))

ross_queue_plot = plot_ross_queue(ross_result.trace; title = "Ross model: repair queue and load")
ross_queue_plot

#jl savefig(ross_queue_plot, plotsdir("ross_queue_and_load.png"))

ross_crash_hist = plot_ross_crash_distribution(ross_runs; title = "Ross model: crash-time distribution")
ross_crash_hist

#jl savefig(ross_crash_hist, plotsdir("ross_crash_histogram.png"))

ross_repairer_plot = plot_ross_comparison(
    ross_repairers;
    xcol = :repairers,
    title = "Ross model: effect of repairer count",
)
ross_repairer_plot

#jl savefig(ross_repairer_plot, plotsdir("ross_repairer_comparison.png"))

# ## Краткие выводы
#
# 1. Для `M/M/c` симуляция хорошо согласуется с аналитическими метриками Эрланга C.
# 2. В модели Росса увеличение числа ремонтников повышает среднее время до отказа.
# 3. Мониторинг длины очереди и загрузки ремонтников позволяет увидеть, когда ремонт становится узким местом системы.
