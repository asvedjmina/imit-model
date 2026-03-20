using DrWatson
@quickactivate "lab_03_models"
using Agents
using DataFrames
using CairoMakie
import StatsBase
using Random

include(srcdir("daisyworld.jl"))

model = daisyworld()

daisycolor(a::Daisy) = a.breed

plotkwargs = (
    agent_color=daisycolor, agent_size = 20, agent_marker = '*',
    heatarray = :temperature,
    heatkwargs = (colorrange = (-20, 60),),
)

plt1, _ = abmplot(model; plotkwargs...)

save(plotsdir("daisy_step001.png"), plt1)

step!(model, 5)
plt2, _ = abmplot(model; heatarray = model.temperature,
    plotkwargs...)

save(plotsdir("daisy_step005.png"), plt2)

step!(model, 40)
plt3, _ = abmplot(model; heatarray = model.temperature,
    plotkwargs...)

save(plotsdir("daisy_step040.png"), plt3)

black(a) = a.breed == :black
white(a) = a.breed == :white
adata = [(black, count), (white, count)]

model2 = daisyworld(; solar_luminosity = 1.0)

agent_df, model_df = run!(model2, 1000; adata)
figure_count = Figure(size = (600, 400));

ax = figure_count[1, 1] = Axis(figure_count, xlabel = "tick", ylabel = "daisy count")
blackl = lines!(ax, agent_df[!, :time], agent_df[!, :count_black],
    color = :black)
whitel = lines!(ax, agent_df[!, :time], agent_df[!, :count_white],
    color = :orange)
Legend(figure_count[1, 2], [blackl, whitel], ["black", "white"], labelsize = 12)
figure_count

save(plotsdir("daisy_count.png"), figure_count)

model3 = daisyworld(solar_luminosity = 1.0, scenario = :ramp)

temperature(model) = StatsBase.mean(model.temperature)
mdata = [temperature, :solar_luminosity]

agent_df3, model_df3 = run!(model3, 1000; adata = adata, mdata = mdata)

figure_lum = CairoMakie.Figure(size = (600, 600));
ax1 = figure_lum[1, 1] = Axis(figure_lum, ylabel = "daisy count")
blackl3 = lines!(ax1, agent_df3[!, :time], agent_df3[!, :count_black],
    color = :red)
whitel3 = lines!(ax1, agent_df3[!, :time], agent_df3[!, :count_white],
    color = :blue)
figure_lum[1, 2] = Legend(figure_lum, [blackl3, whitel3], ["black", "white"])

ax2 = figure_lum[2, 1] = Axis(figure_lum, ylabel = "temperature")
ax3 = figure_lum[3, 1] = Axis(figure_lum, xlabel = "tick", ylabel = "luminosity")
lines!(ax2, model_df3[!, :time], model_df3[!, :temperature], color = :red)
lines!(ax3, model_df3[!, :time], model_df3[!, :solar_luminosity], color = :red)
for ax in (ax1, ax2); ax.xticklabelsvisible = false; end
figure_lum

save(plotsdir("daisy_luminosity.png"), figure_lum)
