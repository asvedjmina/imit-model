# # SIR model: heterogeneity effect
#
# Этот literate-скрипт исследует влияние неоднородных коэффициентов заражения
# по городам. Сравниваются однородный и неоднородный режимы передачи инфекции.

using DrWatson
@quickactivate "lab_04_models"

ENV["GKSwstype"] = "100"

include(srcdir("sir_analysis.jl"))

# ## Сравнение глобальной динамики
#
# Функция `heterogeneity_experiment()` запускает два сценария:
# - однородный `beta`;
# - неоднородный `beta` по трём городам.

base_global, base_city, hetero_global, hetero_city, summary_df, city_summary =
    heterogeneity_experiment()

summary_df

global_plot = plot_heterogeneity_global(base_global, hetero_global)
global_plot

#jl savefig(global_plot, plotsdir("heterogeneity_global_comparison.png"))

# ## Динамика по городам
#
# Здесь нас интересует, как распределяется инфекционная нагрузка между городами.

city_plot = plot_city_dynamics(hetero_city; title = "Неоднородный beta по городам")
city_plot

#jl savefig(city_plot, plotsdir("heterogeneity_city_dynamics.png"))

# ## Составной график

combined_plot = plot(global_plot, city_plot; layout = (2, 1), size = (900, 800))
combined_plot

#jl CSV.write(datadir("heterogeneity_summary.csv"), summary_df)
#jl CSV.write(datadir("heterogeneity_city_summary.csv"), city_summary)
#jl savefig(combined_plot, plotsdir("heterogeneity_effect.png"))

# ## Выводы
#
# 1. Неоднородные значения `beta` меняют распределение нагрузки между городами.
# 2. Из-за миграции различия не изолируются в одном городе.
# 3. Суммарная смертность может даже вырасти, несмотря на неодинаковые параметры.
