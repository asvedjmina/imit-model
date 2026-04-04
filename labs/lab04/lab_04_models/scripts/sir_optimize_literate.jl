# # SIR model: constrained optimization
#
# Этот literate-скрипт решает задачу оптимизации: найти такие параметры,
# которые минимизируют число умерших при ограничении на средний пик
# заболеваемости `mean_peak <= 30%`.

using DrWatson
@quickactivate "lab_04_models"

ENV["GKSwstype"] = "100"

include(srcdir("sir_analysis.jl"))

# ## Поиск допустимых решений
#
# Используется многокритериальная оптимизация `BlackBoxOptim` (`:borg_moea`).
# Оптимизируются параметры:
# - `beta_und`;
# - `detection_time`;
# - `death_rate`.

peak_limit = 0.30
all_df, feasible_df, best = optimize_under_peak_constraint(; peak_limit)

first(feasible_df, min(10, nrow(feasible_df)))

# ## График допустимых решений

opt_plot = plot_constrained_optimization(all_df, feasible_df; peak_limit)
opt_plot

#jl @save datadir("optimization_result.jld2") best feasible_df all_df
#jl CSV.write(datadir("optimization_samples.csv"), all_df)
#jl CSV.write(datadir("optimization_pareto.csv"), feasible_df)
#jl savefig(opt_plot, plotsdir("optimization_pareto.png"))

# ## Выводы
#
# 1. Алгоритм строит Pareto-фронт по двум критериям: средний пик и среднее число умерших.
# 2. Ограничение `mean_peak <= 30%` позволяет выбрать допустимые решения из фронта.
# 3. Лучшие решения смещаются к меньшей заразности, более раннему выявлению и меньшей смертности.
