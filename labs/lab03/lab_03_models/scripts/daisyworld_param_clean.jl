using DrWatson
@quickactivate "lab_03_models"
using Agents
using DataFrames
using CairoMakie
import StatsBase
using Random

include(srcdir("daisyworld.jl"))

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

    save(plotsdir(plt1_name), plt1)
    save(plotsdir(plt2_name), plt2)
    save(plotsdir(plt3_name), plt3)
end

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

    save(plotsdir(plt_name), figure)
end

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
    save(plotsdir(plt_name), figure)
end
