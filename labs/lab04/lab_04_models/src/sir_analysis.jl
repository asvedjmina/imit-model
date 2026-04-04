include("sir_model.jl")

using BlackBoxOptim
using CSV
using DataFrames
using JLD2
using Plots
using Random
using Statistics

function default_parameters()
    return Dict(
        :Ns => [1000, 1000, 1000],
        :beta_und => [0.5, 0.5, 0.5],
        :beta_det => [0.05, 0.05, 0.05],
        :infection_period => 14,
        :detection_time => 7,
        :death_rate => 0.02,
        :reinfection_probability => 0.1,
        :Is => [0, 0, 1],
        :seed => 42,
        :n_steps => 100,
    )
end

function collect_time_series(model, n_steps::Int)
    times = Int[]
    susceptible = Int[]
    infected = Int[]
    recovered = Int[]
    total = Int[]
    city_infected = [Int[] for _ in 1:model.C]

    for step in 1:n_steps
        safe_step!(model)
        push!(times, step)
        push!(susceptible, susceptible_count(model))
        push!(infected, infected_count(model))
        push!(recovered, recovered_count(model))
        push!(total, total_count(model))
        for city in 1:model.C
            push!(city_infected[city], city_status_count(model, city, :I))
        end
    end

    global_df = DataFrame(
        time = times,
        susceptible = susceptible,
        infected = infected,
        recovered = recovered,
        total = total,
    )
    city_df = DataFrame(time = times)
    for city in 1:model.C
        city_df[!, Symbol("city_$(city)_infected")] = city_infected[city]
    end

    return global_df, city_df
end

function run_basic_experiment(params::AbstractDict = default_parameters())
    model = initialize_sir(; pairs(params)...)
    global_df, _ = collect_time_series(model, params[:n_steps])
    agent_df = select(global_df, :time, :susceptible, :infected, :recovered)
    model_df = select(global_df, :time, :total)
    return agent_df, model_df
end

function plot_basic_dynamics(agent_df::DataFrame, model_df::DataFrame)
    plt = plot(
        agent_df.time,
        agent_df.susceptible;
        label = "Восприимчивые",
        xlabel = "Дни",
        ylabel = "Количество",
        linewidth = 2,
    )
    plot!(plt, agent_df.time, agent_df.infected; label = "Инфицированные", linewidth = 2)
    plot!(plt, agent_df.time, agent_df.recovered; label = "Выздоровевшие", linewidth = 2)
    plot!(plt, model_df.time, model_df.total; label = "Всего", linestyle = :dash, linewidth = 2)
    return plt
end

reproduction_number(params::AbstractDict) = first(params[:beta_und]) * params[:infection_period]

function basic_summary(params::AbstractDict = default_parameters())
    model = initialize_sir(; pairs(params)...)
    global_df, _ = collect_time_series(model, params[:n_steps])
    initial_population = sum(params[:Ns])
    peak_idx = argmax(global_df.infected)
    peak_infected = global_df.infected[peak_idx]

    return (
        theoretical_r0 = reproduction_number(params),
        peak_infected = peak_infected,
        peak_fraction = peak_infected / initial_population,
        peak_time = global_df.time[peak_idx],
        final_recovered = global_df.recovered[end],
        deaths = initial_population - global_df.total[end],
    )
end

function beta_metrics(params::AbstractDict)
    model = initialize_sir(; pairs(params)...)
    peak_infected = 0.0
    for _ in 1:params[:n_steps]
        safe_step!(model)
        peak_infected = max(peak_infected, infected_fraction(model; denominator = sum(params[:Ns])))
    end

    initial_population = sum(params[:Ns])
    final_infected = infected_count(model) / initial_population
    final_recovered = recovered_count(model) / initial_population
    total_deaths = initial_population - total_count(model)

    return (
        peak = peak_infected,
        final_inf = final_infected,
        final_rec = final_recovered,
        deaths = total_deaths,
    )
end

function scan_beta_experiment(;
    beta_range = 0.1:0.1:1.0,
    seeds = [42, 43, 44],
    base_params = default_parameters(),
)
    results = NamedTuple[]

    for beta in beta_range
        for seed in seeds
            params = merge(
                copy(base_params),
                Dict(
                    :beta_und => fill(beta, length(base_params[:Ns])),
                    :beta_det => fill(beta / 10, length(base_params[:Ns])),
                    :seed => seed,
                ),
            )
            metrics = beta_metrics(params)
            push!(
                results,
                (
                    beta = beta,
                    seed = seed,
                    peak = metrics.peak,
                    final_inf = metrics.final_inf,
                    final_rec = metrics.final_rec,
                    deaths = metrics.deaths,
                ),
            )
        end
    end

    df = DataFrame(results)
    grouped = combine(
        groupby(df, :beta),
        :peak => mean => :mean_peak,
        :final_inf => mean => :mean_final_inf,
        :final_rec => mean => :mean_final_rec,
        :deaths => mean => :mean_deaths,
    )
    return df, grouped
end

function plot_beta_scan(grouped::DataFrame; population = 3000)
    plt = plot(
        grouped.beta,
        grouped.mean_peak;
        label = "Пик эпидемии",
        xlabel = "Коэффициент заразности beta",
        ylabel = "Доля населения",
        marker = :circle,
        linewidth = 2,
    )
    plot!(plt, grouped.beta, grouped.mean_final_inf; label = "Финальная доля I", marker = :square, linewidth = 2)
    plot!(plt, grouped.beta, grouped.mean_deaths ./ population; label = "Доля умерших", marker = :diamond, linewidth = 2)
    return plt
end

function threshold_study(;
    beta_range = 0.05:0.01:0.30,
    seeds = [42, 43, 44, 45, 46],
    base_params = default_parameters(),
    epidemic_threshold = 0.05,
)
    raw_df, grouped = scan_beta_experiment(; beta_range, seeds, base_params)
    observed_idx = findfirst(grouped.mean_peak .> epidemic_threshold)
    observed_beta = isnothing(observed_idx) ? missing : grouped.beta[observed_idx]
    theoretical_beta = 1 / base_params[:infection_period]

    summary = (
        epidemic_threshold = epidemic_threshold,
        observed_beta = observed_beta,
        theoretical_beta = theoretical_beta,
        threshold_gap = ismissing(observed_beta) ? missing : observed_beta - theoretical_beta,
    )

    return raw_df, grouped, summary
end

function plot_threshold_study(
    grouped::DataFrame;
    epidemic_threshold = 0.05,
    theoretical_beta = 1 / 14,
    observed_beta = missing,
)
    plt = plot(
        grouped.beta,
        grouped.mean_peak;
        xlabel = "beta",
        ylabel = "Средний пик доли I",
        linewidth = 2,
        marker = :circle,
        label = "Средний пик I",
    )
    hline!(plt, [epidemic_threshold]; linestyle = :dash, color = :firebrick, label = "Порог эпидемии 5%")
    vline!(plt, [theoretical_beta]; linestyle = :dashdot, color = :black, label = "Теоретический порог")
    if !ismissing(observed_beta)
        vline!(plt, [observed_beta]; linestyle = :dot, color = :seagreen, label = "Наблюдаемый порог")
    end
    return plt
end

function migration_metrics(params::AbstractDict)
    model = initialize_sir(; pairs(params)...)
    peak_value = 0.0
    peak_time = 0

    for step in 1:params[:n_steps]
        safe_step!(model)
        current_value = infected_fraction(model; denominator = sum(params[:Ns]))
        if current_value > peak_value
            peak_value = current_value
            peak_time = step
        end
    end

    return (peak_time = peak_time, peak_value = peak_value)
end

function scan_migration_experiment(;
    migration_intensities = 0.0:0.1:0.5,
    seeds = [42, 43, 44],
    base_params = default_parameters(),
)
    results = NamedTuple[]
    C = length(base_params[:Ns])

    for intensity in migration_intensities
        for seed in seeds
            params = merge(
                copy(base_params),
                Dict(
                    :migration_rates => create_migration_matrix(C, intensity),
                    :seed => seed,
                    :Is => [1, 0, 0],
                    :n_steps => 150,
                ),
            )
            metrics = migration_metrics(params)
            push!(
                results,
                (
                    migration_intensity = intensity,
                    seed = seed,
                    peak_time = metrics.peak_time,
                    peak_value = metrics.peak_value,
                ),
            )
        end
    end

    df = DataFrame(results)
    grouped = combine(
        groupby(df, :migration_intensity),
        :peak_time => mean => :mean_peak_time,
        :peak_value => mean => :mean_peak_value,
    )
    return df, grouped
end

function plot_migration_effect(grouped::DataFrame; population = 3000)
    p1 = plot(
        grouped.migration_intensity,
        grouped.mean_peak_time;
        marker = :circle,
        linewidth = 2,
        xlabel = "Интенсивность миграции",
        ylabel = "Время до пика, дни",
        label = "Время пика",
    )
    p2 = plot(
        grouped.migration_intensity,
        grouped.mean_peak_value .* population;
        marker = :square,
        linewidth = 2,
        xlabel = "Интенсивность миграции",
        ylabel = "Число инфицированных",
        label = "Пиковая заболеваемость",
    )
    return plot(p1, p2; layout = (2, 1), size = (900, 700))
end

function migration_summary(grouped::DataFrame)
    best_idx = argmin(grouped.mean_peak_time)
    best_row = grouped[best_idx, :]
    return (
        best_intensity = best_row.migration_intensity,
        best_peak_time = best_row.mean_peak_time,
        best_peak_value = best_row.mean_peak_value,
    )
end

function comprehensive_analysis_plot(df::DataFrame)
    grouped = combine(
        groupby(df, :beta),
        :peak => mean => :mean_peak,
        :final_inf => mean => :mean_final_inf,
        :final_rec => mean => :mean_final_rec,
        :deaths => mean => :mean_deaths,
    )

    p1 = plot(
        grouped.beta,
        grouped.mean_peak;
        label = "Пик",
        xlabel = "beta",
        ylabel = "Доля инфицированных",
        linewidth = 2,
    )
    plot!(p1, grouped.beta, grouped.mean_final_inf; label = "Финальная доля I", linewidth = 2)

    p2 = plot(
        grouped.beta,
        grouped.mean_deaths;
        label = "Умершие",
        xlabel = "beta",
        ylabel = "Число умерших",
        linewidth = 2,
        color = :firebrick,
    )

    p3 = plot(
        grouped.beta,
        grouped.mean_final_rec;
        label = "Доля выздоровевших",
        xlabel = "beta",
        ylabel = "Доля",
        linewidth = 2,
        color = :seagreen,
    )

    return plot(p1, p2, p3; layout = (3, 1), size = (900, 1000))
end

function scenario_parameters()
    return [
        Dict(
            :label => "baseline",
            :beta_und => [0.5, 0.5, 0.5],
            :beta_det => [0.05, 0.05, 0.05],
            :migration_intensity => 0.1,
            :detection_time => 7,
            :n_steps => 120,
            :seed => 42,
        ),
        Dict(
            :label => "heterogeneous_beta",
            :beta_und => [0.75, 0.45, 0.25],
            :beta_det => [0.08, 0.05, 0.03],
            :migration_intensity => 0.1,
            :detection_time => 7,
            :n_steps => 120,
            :seed => 52,
        ),
        Dict(
            :label => "fast_migration",
            :beta_und => [0.5, 0.5, 0.5],
            :beta_det => [0.05, 0.05, 0.05],
            :migration_intensity => 0.4,
            :detection_time => 7,
            :n_steps => 120,
            :seed => 62,
        ),
        Dict(
            :label => "slow_detection",
            :beta_und => [0.5, 0.5, 0.5],
            :beta_det => [0.05, 0.05, 0.05],
            :migration_intensity => 0.1,
            :detection_time => 10,
            :n_steps => 120,
            :seed => 72,
        ),
    ]
end

function run_scenario(spec::AbstractDict; base_params = default_parameters())
    params = merge(
        copy(base_params),
        Dict(
            :beta_und => spec[:beta_und],
            :beta_det => spec[:beta_det],
            :detection_time => spec[:detection_time],
            :seed => spec[:seed],
            :n_steps => spec[:n_steps],
            :migration_rates => create_migration_matrix(length(base_params[:Ns]), spec[:migration_intensity]),
        ),
    )

    model = initialize_sir(; pairs(params)...)
    global_df, city_df = collect_time_series(model, params[:n_steps])
    initial_population = sum(params[:Ns])
    peak_time = global_df.time[argmax(global_df.infected)]

    summary = (
        scenario = spec[:label],
        peak_infected = maximum(global_df.infected),
        peak_fraction = maximum(global_df.infected) / initial_population,
        peak_time = peak_time,
        final_recovered = global_df.recovered[end],
        deaths = initial_population - global_df.total[end],
    )

    return global_df, city_df, summary
end

function plot_scenario(global_df::DataFrame, city_df::DataFrame; title = "")
    p1 = plot(
        global_df.time,
        global_df.susceptible;
        label = "S",
        linewidth = 2,
        xlabel = "Дни",
        ylabel = "Количество",
        title = title,
    )
    plot!(p1, global_df.time, global_df.infected; label = "I", linewidth = 2)
    plot!(p1, global_df.time, global_df.recovered; label = "R", linewidth = 2)

    p2 = plot(
        city_df.time,
        city_df.city_1_infected;
        label = "Город 1",
        linewidth = 2,
        xlabel = "Дни",
        ylabel = "Инфицированные",
    )
    plot!(p2, city_df.time, city_df.city_2_infected; label = "Город 2", linewidth = 2)
    plot!(p2, city_df.time, city_df.city_3_infected; label = "Город 3", linewidth = 2)

    return plot(p1, p2; layout = (2, 1), size = (900, 700))
end

function heterogeneity_experiment(;
    base_params = default_parameters(),
    n_steps = 120,
    seed = 52,
)
    baseline_spec = Dict(
        :label => "homogeneous_beta",
        :beta_und => [0.5, 0.5, 0.5],
        :beta_det => [0.05, 0.05, 0.05],
        :migration_intensity => 0.1,
        :detection_time => 7,
        :n_steps => n_steps,
        :seed => seed,
    )
    hetero_spec = Dict(
        :label => "heterogeneous_beta",
        :beta_und => [0.80, 0.45, 0.20],
        :beta_det => [0.08, 0.05, 0.02],
        :migration_intensity => 0.1,
        :detection_time => 7,
        :n_steps => n_steps,
        :seed => seed,
    )

    base_global, base_city, base_summary = run_scenario(baseline_spec; base_params)
    hetero_global, hetero_city, hetero_summary = run_scenario(hetero_spec; base_params)

    city_summary = DataFrame(
        city = 1:3,
        beta_und = hetero_spec[:beta_und],
        peak_infected = [
            maximum(hetero_city[!, :city_1_infected]),
            maximum(hetero_city[!, :city_2_infected]),
            maximum(hetero_city[!, :city_3_infected]),
        ],
    )
    summary_df = DataFrame([base_summary, hetero_summary])

    return base_global, base_city, hetero_global, hetero_city, summary_df, city_summary
end

function plot_heterogeneity_global(base_global::DataFrame, hetero_global::DataFrame)
    plt = plot(
        base_global.time,
        base_global.infected;
        label = "Однородный beta",
        xlabel = "Дни",
        ylabel = "Число инфицированных",
        linewidth = 2,
    )
    plot!(plt, hetero_global.time, hetero_global.infected; label = "Неоднородный beta", linewidth = 2)
    return plt
end

function plot_city_dynamics(city_df::DataFrame; title = "Динамика по городам")
    plt = plot(
        city_df.time,
        city_df.city_1_infected;
        label = "Город 1",
        xlabel = "Дни",
        ylabel = "Инфицированные",
        linewidth = 2,
        title = title,
    )
    plot!(plt, city_df.time, city_df.city_2_infected; label = "Город 2", linewidth = 2)
    plot!(plt, city_df.time, city_df.city_3_infected; label = "Город 3", linewidth = 2)
    return plt
end

function quarantine_experiment(;
    thresholds = [0.05, 0.10, 0.15, 0.20],
    base_params = default_parameters(),
    beta = 0.20,
    migration_intensity = 0.10,
    n_steps = 120,
    seed = 91,
)
    results = NamedTuple[]
    dynamics = Dict{String, NamedTuple}()
    C = length(base_params[:Ns])

    scenarios = vcat([nothing], thresholds)
    for threshold in scenarios
        params = merge(
            copy(base_params),
            Dict(
                :beta_und => fill(beta, C),
                :beta_det => fill(beta / 10, C),
                :migration_rates => create_migration_matrix(C, migration_intensity),
                :Is => [1, 0, 0],
                :seed => seed,
                :n_steps => n_steps,
                :quarantine_threshold => threshold,
            ),
        )

        model = initialize_sir(; pairs(params)...)
        global_df, city_df = collect_time_series(model, params[:n_steps])
        initial_population = sum(params[:Ns])
        peak_idx = argmax(global_df.infected)
        label = isnothing(threshold) ? "no_quarantine" : "q_" * string(round(Int, threshold * 100))

        push!(
            results,
            (
                scenario = label,
                quarantine_threshold = something(threshold, -1.0),
                peak_fraction = global_df.infected[peak_idx] / initial_population,
                peak_time = global_df.time[peak_idx],
                deaths = initial_population - global_df.total[end],
                locked_cities = locked_city_count(model),
            ),
        )
        dynamics[label] = (global_df = global_df, city_df = city_df)
    end

    summary_df = DataFrame(results)
    return summary_df, dynamics
end

function plot_quarantine_summary(summary_df::DataFrame)
    labels = replace.(summary_df.scenario, "q_" => "q=")
    p1 = bar(
        labels,
        summary_df.peak_fraction;
        label = "Пиковая доля I",
        xlabel = "Сценарий",
        ylabel = "Доля",
    )
    hline!(p1, [0.30]; linestyle = :dash, color = :black, label = "Ограничение 30%")

    p2 = bar(
        labels,
        summary_df.deaths;
        label = "Смерти",
        xlabel = "Сценарий",
        ylabel = "Количество",
        color = :firebrick,
    )

    return plot(p1, p2; layout = (2, 1), size = (900, 700))
end

function plot_quarantine_dynamics(baseline_df::DataFrame, quarantine_df::DataFrame; quarantine_label = "")
    plt = plot(
        baseline_df.time,
        baseline_df.infected;
        label = "Без карантина",
        xlabel = "Дни",
        ylabel = "Инфицированные",
        linewidth = 2,
    )
    plot!(plt, quarantine_df.time, quarantine_df.infected; label = quarantine_label, linewidth = 2)
    return plt
end

function evaluate_policy(beta_und, detection_time, death_rate;
    base_params = default_parameters(),
    replicates = 5,
    n_steps = 100,
)
    peak_values = Float64[]
    death_values = Float64[]

    for rep in 1:replicates
        params = merge(
            copy(base_params),
            Dict(
                :beta_und => fill(beta_und, length(base_params[:Ns])),
                :beta_det => fill(beta_und / 10, length(base_params[:Ns])),
                :detection_time => detection_time,
                :death_rate => death_rate,
                :seed => 500 + rep,
                :n_steps => n_steps,
            ),
        )
        metrics = beta_metrics(params)
        push!(peak_values, metrics.peak)
        push!(death_values, metrics.deaths)
    end

    return mean(peak_values), mean(death_values)
end

function optimization_cost(x;
    base_params = default_parameters(),
    replicates = 5,
    n_steps = 100,
)
    beta_und = clamp(x[1], 0.1, 1.0)
    detection_time = clamp(round(Int, x[2]), 3, 14)
    death_rate = clamp(x[3], 0.01, 0.1)

    mean_peak, mean_deaths = evaluate_policy(
        beta_und,
        detection_time,
        death_rate;
        base_params,
        replicates,
        n_steps,
    )

    return (mean_peak, mean_deaths)
end

function pareto_frontier_dataframe(result, peak_limit)
    rows = NamedTuple[]

    for individual in pareto_frontier(result)
        candidate = params(individual)
        objectives = archived_fitness(individual)
        objective_tuple = objectives.orig
        detection_time = clamp(round(Int, candidate[2]), 3, 14)

        push!(
            rows,
            (
                beta_und = candidate[1],
                detection_time = detection_time,
                death_rate = candidate[3],
                mean_peak = objective_tuple[1],
                mean_deaths = objective_tuple[2],
                feasible = objective_tuple[1] <= peak_limit,
            ),
        )
    end

    frontier_df = DataFrame(rows)
    nrow(frontier_df) == 0 && return frontier_df

    sort!(frontier_df, [:mean_peak, :mean_deaths, :beta_und, :detection_time, :death_rate])
    unique!(frontier_df)

    return frontier_df
end

function optimize_under_peak_constraint(;
    peak_limit = 0.30,
    max_evals = 30,
    replicates = 2,
    n_steps = 60,
    seed = 2026,
    base_params = default_parameters(),
)
    Random.seed!(seed)

    result = bboptimize(
        x -> optimization_cost(
            x;
            base_params,
            replicates,
            n_steps,
        ),
        Method = :borg_moea,
        FitnessScheme = ParetoFitnessScheme{2}(is_minimizing = true),
        SearchRange = [(0.1, 1.0), (3.0, 14.0), (0.01, 0.1)],
        NumDimensions = 3,
        MaxFuncEvals = max_evals,
        TraceMode = :compact,
    )

    all_df = pareto_frontier_dataframe(result, peak_limit)
    feasible_df = filter(:feasible => identity, all_df)
    sort!(feasible_df, [:mean_deaths, :mean_peak, :beta_und, :detection_time, :death_rate])

    best = nrow(feasible_df) == 0 ? nothing : feasible_df[1, :]
    return all_df, feasible_df, best
end

function plot_constrained_optimization(all_df::DataFrame, feasible_df::DataFrame; peak_limit = 0.30)
    plt = scatter(
        all_df.mean_peak,
        all_df.mean_deaths;
        xlabel = "Средний пик доли I",
        ylabel = "Среднее число умерших",
        label = "Точки Парето-фронта",
        alpha = 0.65,
        markerstrokewidth = 0,
        color = :steelblue,
    )
    if nrow(feasible_df) > 0
        scatter!(plt, feasible_df.mean_peak, feasible_df.mean_deaths; label = "Допустимые", color = :darkgreen)
    end
    vline!(plt, [peak_limit]; linestyle = :dash, color = :firebrick, label = "Предел 30%")
    return plt
end

function scenario_comparison_plot(summary_df::DataFrame)
    p1 = bar(
        summary_df.scenario,
        summary_df.peak_fraction;
        label = "Пиковая доля I",
        xlabel = "Сценарий",
        ylabel = "Доля",
        legend = :topright,
    )
    p2 = bar(
        summary_df.scenario,
        summary_df.deaths;
        label = "Число умерших",
        xlabel = "Сценарий",
        ylabel = "Количество",
        color = :firebrick,
        legend = :topright,
    )
    return plot(p1, p2; layout = (2, 1), size = (900, 700))
end
