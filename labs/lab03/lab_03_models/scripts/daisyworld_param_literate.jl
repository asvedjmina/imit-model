# # Модель Daisyworld: исследование параметров
#
# Данный скрипт расширяет базовую модель Daisyworld, исследуя влияние
# различных параметров модели на её поведение. Мы варьируем максимальный
# возраст маргариток и начальную долю белых маргариток.
#
# ## Параметрическое пространство
#
# | Параметр | Значения | Описание |
# |----------|----------|---------|
# | `max_age` | 25, 40 | Максимальный возраст маргаритки |
# | `init_white` | 0.2, 0.8 | Начальная доля белых маргариток |
# | `init_black` | 0.2 | Начальная доля чёрных маргариток |
#
# Всего получается 4 комбинации параметров (2×2).

# ## Подключение пакетов

using DrWatson
@quickactivate "lab_03_models"
using Agents
using DataFrames
using CairoMakie
import StatsBase
using Random

# ## Загрузка модели

include(srcdir("daisyworld.jl"))

# ## Базовая визуализация с параметрами
#
# Для каждой комбинации параметров создаём модель и визуализируем
# состояние на шагах 0, 5 и 45.

param_dict = Dict(
    :griddims => (30, 30),
    :max_age => [25, 40],
    :init_white => [0.2, 0.8],
    :init_black => 0.2,
    :albedo_white => 0.75,
    :albedo_black => 0.25,
    :surface_albedo => 0.4,
    :solar_change => 0.005,
    :solar_luminosity => 1.0,
    :scenario => :default,
    :seed => 165,
)

params_list = dict_list(param_dict)

for params in params_list

    model = daisyworld(;params...)

    daisycolor(a::Daisy) = a.breed

    plotkwargs = (
        agent_color=daisycolor,
        agent_size = 20,
        agent_marker = '*',
        heatarray = :temperature,
        heatkwargs = (colorrange = (-20, 60),),
    )
    plt1, _ = abmplot(model; plotkwargs...)

    step!(model, 5)
    plt2, _ = abmplot(model; heatarray = model.temperature,
        plotkwargs...)

    step!(model, 40)
    plt3, _ = abmplot(model; heatarray = model.temperature,
        plotkwargs...)

    plt1_name = savename("daisyworld",params) * "_step01" * ".png"
    plt2_name = savename("daisyworld",params) * "_step04" * ".png"
    plt3_name = savename("daisyworld",params) * "_step40" * ".png"

#jl     save(plotsdir(plt1_name), plt1)
#jl     save(plotsdir(plt2_name), plt2)
#jl     save(plotsdir(plt3_name), plt3)
end

# ## Динамика числа маргариток с параметрами
#
# Построим графики изменения числа маргариток для каждой комбинации параметров.

black(a) = a.breed == :black
white(a) = a.breed == :white
adata = [(black, count), (white, count)]

for params in params_list

    model = daisyworld(;params...)

    agent_df, model_df = run!(model, 1000; adata)
    figure = Figure(size = (600, 400));
    ax = figure[1, 1] = Axis(figure, xlabel = "tick", ylabel = "daisy count")
    blackl = lines!(ax, agent_df[!, :time], agent_df[!, :count_black],
        color = :black)
    whitel = lines!(ax, agent_df[!, :time], agent_df[!, :count_white],
        color = :orange)
    Legend(figure[1, 2], [blackl, whitel], ["black", "white"],
        labelsize = 12)

    plt_name = savename("daisy-count",params) * ".png"

#jl     save(plotsdir(plt_name), figure)
end

# ## Динамика модели с параметрами (ramp-сценарий)
#
# Исследуем влияние параметров на поведение модели при изменении светимости.

param_dict_ramp = Dict(
    :griddims => (30, 30),
    :max_age => [25, 40],
    :init_white => [0.2, 0.8],
    :init_black => 0.2,
    :albedo_white => 0.75,
    :albedo_black => 0.25,
    :surface_albedo => 0.4,
    :solar_change => 0.005,
    :solar_luminosity => 1.0,
    :scenario => :ramp,
    :seed => 165,
)

params_list_ramp = dict_list(param_dict_ramp)

for params in params_list_ramp

    model = daisyworld(;params...)

    temperature(model) = StatsBase.mean(model.temperature)
    mdata = [temperature, :solar_luminosity]

    agent_df, model_df = run!(model, 1000; adata = adata, mdata = mdata)

    figure = CairoMakie.Figure(size = (600, 600));
    ax1 = figure[1, 1] = Axis(figure, ylabel = "daisy count")
    blackl = lines!(ax1, agent_df[!, :time], agent_df[!, :count_black],
        color = :red)
    whitel = lines!(ax1, agent_df[!, :time], agent_df[!, :count_white],
        color = :blue)
    figure[1, 2] = Legend(figure, [blackl, whitel], ["black", "white"])

    ax2 = figure[2, 1] = Axis(figure, ylabel = "temperature")
    ax3 = figure[3, 1] = Axis(figure, xlabel = "tick", ylabel = "luminosity")
    lines!(ax2, model_df[!, :time], model_df[!, :temperature], color = :red)
    lines!(ax3, model_df[!, :time], model_df[!, :solar_luminosity], color = :red)
    for ax in (ax1, ax2); ax.xticklabelsvisible = false; end

    plt_name = savename("daisy-luminosity",params) * ".png"
#jl     save(plotsdir(plt_name), figure)
end

# ## Выводы по параметрическому исследованию
#
# 1. **Максимальный возраст (`max_age`)**: увеличение возраста маргариток
#    делает популяцию более устойчивой, так как маргаритки живут дольше
#    и успевают произвести больше потомства.
# 2. **Начальная доля белых (`init_white`)**: при высокой начальной доле
#    белых маргариток (0.8) температура планеты изначально ниже,
#    что влияет на баланс популяций.
# 3. **Взаимодействие параметров**: комбинация высокого возраста и
#    большой начальной доли белых создаёт наиболее устойчивую систему.
