# # Модель SIR в терминах сети Петри
#
# В этой лабораторной работе эпидемическая модель SIR
# рассматривается как размеченная сеть Петри.
# Это позволяет единым образом описать структуру переходов,
# стохастическую динамику и детерминированную аппроксимацию.
#
# ## Идея модели
#
# В системе есть три позиции:
#
# - `S` — восприимчивые;
# - `I` — инфицированные;
# - `R` — выздоровевшие.
#
# И два перехода:
#
# - `infection`: `S + I -> I + I`;
# - `recovery`: `I -> R`.
#
# Стохастический запуск использует алгоритм Гиллеспи,
# а детерминированный — интегрирование системы SIR
# методом Рунге-Кутты четвёртого порядка на фиксированной сетке.

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

# ## Базовые параметры
#
# Выбирается типичная закрытая популяция из 1000 человек
# с десятью начальными инфицированными.

beta = 0.30
gamma = 0.10
tmax = 100.0
saveat = 0.5
seed = 123

net, u0, states = build_sir_network(beta, gamma)

# ## Схема сети Петри
#
# Для отчёта сохраняются как графическое представление,
# так и текстовое описание в формате Graphviz DOT.

network_plot = plot_sir_network(net)
network_plot

#jl savefig(network_plot, plotsdir("sir_network.png"))
#jl write(joinpath(plotsdir(), "sir_network.dot"), to_graphviz_sir(net))

# ## Детерминированная симуляция
#
# Решается непрерывная версия модели SIR
# на равномерной временной сетке.

df_det = simulate_deterministic(
    net,
    u0,
    (0.0, tmax);
    saveat = saveat,
    rates = [beta, gamma],
)

det_summary = DataFrame([summarize_trajectory(df_det)])
det_summary

#jl CSV.write(datadir("sir_det.csv"), df_det)
#jl CSV.write(datadir("sir_det_summary.csv"), det_summary)

det_plot = plot_sir(df_det; title = "Deterministic SIR dynamics")
det_plot

#jl savefig(det_plot, plotsdir("sir_det_dynamics.png"))

# ## Стохастическая симуляция
#
# Запускается прямой метод Гиллеспи.
# При большой популяции его траектория близка к сглаженной ОДУ-модели,
# но сохраняет случайные флуктуации.

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

#jl CSV.write(datadir("sir_stoch.csv"), df_stoch)
#jl CSV.write(datadir("sir_stoch_summary.csv"), stoch_summary)

stoch_plot = plot_sir(df_stoch; title = "Stochastic SIR dynamics")
stoch_plot

#jl savefig(stoch_plot, plotsdir("sir_stoch_dynamics.png"))

# ## Сравнение по числу инфицированных
#
# В итоговом отчёте полезно сравнивать именно траекторию `I(t)`.

comparison_plot = plot_comparison(df_det, df_stoch)
comparison_plot

#jl savefig(comparison_plot, plotsdir("comparison.png"))

# ## Чувствительность по параметру beta
#
# Исследуется влияние скорости заражения на пик эпидемии
# и конечное число выздоровевших.

beta_range = 0.10:0.05:0.80
df_scan = parameter_scan(beta_range; gamma = gamma, tmax = tmax, saveat = saveat)
df_scan

#jl CSV.write(datadir("sir_scan.csv"), df_scan)

scan_plot = plot_scan(df_scan)
scan_plot

#jl savefig(scan_plot, plotsdir("sir_scan.png"))

sensitivity_plot = plot_sensitivity(df_scan)
sensitivity_plot

#jl savefig(sensitivity_plot, plotsdir("sensitivity.png"))

# ## Анимация детерминированной динамики
#
# GIF отражает, как популяция последовательно
# перераспределяется между тремя компартментами.

df_anim = simulate_deterministic(
    net,
    u0,
    (0.0, tmax);
    saveat = 1.0,
    rates = [beta, gamma],
)

#jl animate_sir(df_anim, plotsdir("sir_animation.gif"); fps = 8)

# ## Краткие выводы
#
# 1. Сеть Петри корректно воспроизводит классическую структуру SIR-модели.
# 2. Стохастическая траектория колеблется вокруг детерминированной,
#    но сохраняет те же основные фазы эпидемии.
# 3. При росте `beta` пик эпидемии и итоговое число переболевших увеличиваются.
