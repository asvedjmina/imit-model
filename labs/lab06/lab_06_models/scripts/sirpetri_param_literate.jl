# # Модель SIR в терминах сети Петри: набор параметров
#
# В этом сценарии выполняется серия детерминированных прогонов
# для набора пар `(beta, gamma)`.
# Это соответствует дополнительной части лабораторной работы,
# где literate-код расширяется расчётом для нескольких конфигураций.

using DrWatson
@quickactivate "lab_06_models"

ENV["GKSwstype"] = "100"

include(joinpath(@__DIR__, "..", "src", "SIRPetri.jl"))
using .SIRPetri
using CSV
using DataFrames
using Plots

mkpath(datadir())
mkpath(plotsdir())

# ## Набор параметров
#
# Значения подобраны так, чтобы охватить
# как слабые, так и выраженные эпидемические режимы.

parameter_pairs = [
    (beta = 0.15, gamma = 0.15),
    (beta = 0.20, gamma = 0.10),
    (beta = 0.25, gamma = 0.10),
    (beta = 0.30, gamma = 0.10),
    (beta = 0.35, gamma = 0.10),
    (beta = 0.40, gamma = 0.08),
    (beta = 0.45, gamma = 0.08),
    (beta = 0.50, gamma = 0.12),
    (beta = 0.55, gamma = 0.12),
]

grid_df = parameter_grid(parameter_pairs; tmax = 100.0, saveat = 0.5)
grid_df

#jl CSV.write(datadir("sir_parameter_grid.csv"), grid_df)

# ## Тепловая карта максимума инфицированных
#
# Такой вид представления быстро показывает,
# в каких областях параметров эпидемия наиболее тяжёлая.

peak_heatmap = plot_parameter_heatmap(grid_df; metric = :peak_I)
peak_heatmap

#jl savefig(peak_heatmap, plotsdir("sir_parameter_heatmap.png"))

# ## Тепловая карта конечного числа выздоровевших

recovered_heatmap = plot_parameter_heatmap(grid_df; metric = :final_R)
recovered_heatmap

#jl savefig(recovered_heatmap, plotsdir("sir_finalR_heatmap.png"))

# ## Сводный ранжированный график
#
# Для отчёта дополнительно удобно отсортировать конфигурации
# по высоте эпидемического пика.

sorted_df = sort(grid_df, :peak_I, rev = true)
labels = ["b=$(row.beta), g=$(row.gamma)" for row in eachrow(sorted_df)]

rank_plot = bar(
    labels,
    sorted_df.peak_I;
    xlabel = "Parameter pair",
    ylabel = "Peak I",
    title = "Peak infected for parameter pairs",
    legend = false,
    xrotation = 35,
    size = (900, 500),
)
rank_plot

#jl savefig(rank_plot, plotsdir("sir_parameter_rank.png"))

# ## Краткие выводы
#
# 1. Отношение `beta / gamma` определяет тяжесть эпидемии.
# 2. При повышении `beta` и снижении `gamma`
#    максимум инфицированных заметно возрастает.
# 3. Набор параметров удобно анализировать через тепловые карты
#    и ранжированные столбчатые диаграммы.
