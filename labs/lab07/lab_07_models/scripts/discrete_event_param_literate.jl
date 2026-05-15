# # Дискретно-событийные модели: параметрические прогоны
#
# В этом сценарии выполняются серии запусков для набора параметров:
#
# - число каналов в `M/M/c`;
# - число машин в модели Росса;
# - число ремонтников в модели Росса.

using DrWatson
@quickactivate "lab_07_models"

ENV["GKSwstype"] = "100"

include(joinpath(@__DIR__, "..", "src", "DiscreteEventModels.jl"))
using .DiscreteEventModels
using CSV
using Plots

mkpath(datadir())
mkpath(plotsdir())

# ## Сканирование по числу каналов в M/M/c
#
# При увеличении `c` вероятность ожидания и средняя задержка
# должны уменьшаться.

mmc_sweep = mmc_server_sweep(1:5; lambda = 0.45, mu = 0.5, num_customers = 1000, seed = 900)
mmc_sweep

#jl CSV.write(datadir("mmc_server_sweep.csv"), mmc_sweep)

mmc_sweep_plot = plot_mmc_sweep(mmc_sweep)
mmc_sweep_plot

#jl savefig(mmc_sweep_plot, plotsdir("mmc_server_sweep.png"))

# ## Сканирование по числу основных машин в модели Росса
#
# Чем больше одновременно работающих машин при фиксированном резерве,
# тем быстрее система исчерпывает запас и тем меньше время до сбоя.

ross_machine_df = ross_machine_sweep(6:2:14; spares = 3, repairers = 2, runs = 35, seed = 1200)
ross_machine_df

#jl CSV.write(datadir("ross_machine_sweep.csv"), ross_machine_df)

ross_machine_plot = plot_ross_comparison(
    ross_machine_df;
    xcol = :machines,
    title = "Ross model: crash time vs number of machines",
)
ross_machine_plot

#jl savefig(ross_machine_plot, plotsdir("ross_machine_sweep.png"))

ross_machine_load = plot_ross_utilization_sweep(
    ross_machine_df;
    xcol = :machines,
    title = "Ross model: repair load vs number of machines",
)
ross_machine_load

#jl savefig(ross_machine_load, plotsdir("ross_machine_load.png"))

# ## Сканирование по числу ремонтников
#
# Дополнительный прогон показывает, насколько
# чувствительна система к расширению ремонтной службы.

ross_repairer_df = ross_repairer_sweep(1:4; machines = 10, spares = 3, runs = 35, seed = 1600)
ross_repairer_df

#jl CSV.write(datadir("ross_repairer_sweep_extended.csv"), ross_repairer_df)

ross_repairer_plot = plot_ross_comparison(
    ross_repairer_df;
    xcol = :repairers,
    title = "Ross model: crash time vs number of repairers",
)
ross_repairer_plot

#jl savefig(ross_repairer_plot, plotsdir("ross_repairer_sweep_extended.png"))

ross_repairer_load = plot_ross_utilization_sweep(
    ross_repairer_df;
    xcol = :repairers,
    title = "Ross model: repair load vs number of repairers",
)
ross_repairer_load

#jl savefig(ross_repairer_load, plotsdir("ross_repairer_load.png"))

# ## Краткие выводы
#
# 1. В `M/M/c` увеличение числа каналов резко снижает вероятность ожидания.
# 2. В модели Росса рост числа работающих машин ускоряет падение системы.
# 3. Добавление ремонтников увеличивает среднее время до отказа и уменьшает очередь на ремонт.
