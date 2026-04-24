module SIRPetri

using DataFrames
using Plots
using Random

export LabelledPetriNet,
    animate_sir,
    build_sir_network,
    parameter_grid,
    parameter_scan,
    plot_comparison,
    plot_parameter_heatmap,
    plot_scan,
    plot_sensitivity,
    plot_sir,
    plot_sir_network,
    sample_path,
    simulate_deterministic,
    simulate_stochastic,
    sir_ode,
    summarize_trajectory,
    to_graphviz_sir

struct LabelledPetriNet
    states::Vector{Symbol}
    transitions::Vector{Symbol}
    inputs::Vector{Vector{Int}}
    outputs::Vector{Vector{Int}}
    rates::Dict{Symbol, Float64}
end

default_rates(net::LabelledPetriNet) = [
    net.rates[:infection],
    net.rates[:recovery],
]

function build_sir_network(beta::Real = 0.3, gamma::Real = 0.1; initial_state = [990.0, 10.0, 0.0])
    states = [:S, :I, :R]
    transitions = [:infection, :recovery]
    inputs = [[1, 2], [2]]
    outputs = [[2, 2], [3]]
    net = LabelledPetriNet(
        states,
        transitions,
        inputs,
        outputs,
        Dict(:infection => float(beta), :recovery => float(gamma)),
    )
    u0 = Float64.(initial_state)
    return net, u0, copy(states)
end

function sir_ode(net::LabelledPetriNet, rates = default_rates(net))
    function f!(du, u, p, t)
        S, I, R = u
        beta, gamma = rates
        population = max(sum(u), 1.0)
        infection_rate = beta * S * I / population
        recovery_rate = gamma * I

        du[1] = -infection_rate
        du[2] = infection_rate - recovery_rate
        du[3] = recovery_rate
        return nothing
    end
    return f!
end

function rk4_step(f!, u::Vector{Float64}, t::Float64, dt::Float64)
    k1 = similar(u)
    k2 = similar(u)
    k3 = similar(u)
    k4 = similar(u)

    f!(k1, u, nothing, t)
    f!(k2, u .+ 0.5 .* dt .* k1, nothing, t + 0.5 * dt)
    f!(k3, u .+ 0.5 .* dt .* k2, nothing, t + 0.5 * dt)
    f!(k4, u .+ dt .* k3, nothing, t + dt)

    next_u = u .+ dt ./ 6.0 .* (k1 .+ 2.0 .* k2 .+ 2.0 .* k3 .+ k4)
    next_u = max.(next_u, 0.0)

    target_total = sum(u)
    current_total = sum(next_u)
    if current_total > 0
        next_u .*= target_total / current_total
    end

    return next_u
end

function time_grid(tspan::Tuple{<:Real, <:Real}, saveat::Real)
    t0, t1 = float(tspan[1]), float(tspan[2])
    times = collect(t0:float(saveat):t1)
    if isempty(times) || times[end] < t1
        push!(times, t1)
    end
    return times
end

function marking_dataframe(times, states)
    return DataFrame(
        time = Float64.(times),
        S = [state[1] for state in states],
        I = [state[2] for state in states],
        R = [state[3] for state in states],
    )
end

function simulate_deterministic(
    net::LabelledPetriNet,
    u0::AbstractVector{<:Real},
    tspan::Tuple{<:Real, <:Real};
    saveat::Real = 0.1,
    rates = default_rates(net),
)
    f! = sir_ode(net, rates)
    times = time_grid(tspan, saveat)
    state = Float64.(u0)
    states = [copy(state)]

    for idx in 1:length(times)-1
        dt = times[idx + 1] - times[idx]
        state = rk4_step(f!, state, times[idx], dt)
        push!(states, copy(state))
    end

    return marking_dataframe(times, states)
end

function simulate_stochastic(
    net::LabelledPetriNet,
    u0::AbstractVector{<:Real},
    tspan::Tuple{<:Real, <:Real};
    rates = default_rates(net),
    rng::AbstractRNG = Random.default_rng(),
)
    beta, gamma = rates
    t0, t1 = float(tspan[1]), float(tspan[2])
    state = Int.(round.(u0))
    t = t0
    times = [t]
    states = [copy(state)]

    while t < t1
        S, I, R = state
        population = max(sum(state), 1)
        infection_rate = beta * S * I / population
        recovery_rate = gamma * I
        total_rate = infection_rate + recovery_rate

        if total_rate <= 0
            break
        end

        dt = -log(rand(rng)) / total_rate
        next_t = t + dt
        if next_t > t1
            break
        end

        if rand(rng) * total_rate < infection_rate
            state[1] -= 1
            state[2] += 1
        else
            state[2] -= 1
            state[3] += 1
        end

        t = next_t
        push!(times, t)
        push!(states, copy(state))
    end

    if times[end] < t1
        push!(times, t1)
        push!(states, copy(state))
    end

    return marking_dataframe(times, states)
end

function summarize_trajectory(df::DataFrame)
    peak_index = argmax(df.I)
    return (
        peak_I = maximum(df.I),
        peak_time = df.time[peak_index],
        final_S = df.S[end],
        final_R = df.R[end],
    )
end

function parameter_scan(
    beta_values;
    gamma::Real = 0.1,
    tmax::Real = 100.0,
    saveat::Real = 0.5,
    initial_state = [990.0, 10.0, 0.0],
)
    rows = NamedTuple[]
    for beta in beta_values
        net, u0, _ = build_sir_network(beta, gamma; initial_state = initial_state)
        df = simulate_deterministic(net, u0, (0.0, tmax); saveat = saveat, rates = [beta, gamma])
        summary = summarize_trajectory(df)
        push!(
            rows,
            (
                beta = float(beta),
                gamma = float(gamma),
                peak_I = summary.peak_I,
                peak_time = summary.peak_time,
                final_R = summary.final_R,
            ),
        )
    end
    return DataFrame(rows)
end

function parameter_grid(
    parameter_pairs;
    tmax::Real = 100.0,
    saveat::Real = 0.5,
    initial_state = [990.0, 10.0, 0.0],
)
    rows = NamedTuple[]
    for pair in parameter_pairs
        beta = pair.beta
        gamma = pair.gamma
        net, u0, _ = build_sir_network(beta, gamma; initial_state = initial_state)
        df = simulate_deterministic(net, u0, (0.0, tmax); saveat = saveat, rates = [beta, gamma])
        summary = summarize_trajectory(df)
        push!(
            rows,
            (
                beta = float(beta),
                gamma = float(gamma),
                peak_I = summary.peak_I,
                peak_time = summary.peak_time,
                final_R = summary.final_R,
                final_S = summary.final_S,
            ),
        )
    end
    return DataFrame(rows)
end

function plot_sir(df::DataFrame; title::AbstractString = "SIR dynamics")
    return plot(
        df.time,
        hcat(df.S, df.I, df.R);
        label = ["S (Susceptible)" "I (Infected)" "R (Recovered)"],
        xlabel = "Time",
        ylabel = "Population",
        linewidth = 2,
        title = title,
    )
end

function plot_scan(df::DataFrame)
    return plot(
        df.beta,
        hcat(df.peak_I, df.final_R);
        label = ["Peak I" "Final R"],
        marker = :circle,
        xlabel = "beta",
        ylabel = "Population",
        linewidth = 2,
        title = "Sensitivity to beta",
    )
end

plot_sensitivity(df::DataFrame) = plot(
    df.beta,
    df.peak_I;
    marker = :circle,
    xlabel = "beta",
    ylabel = "Peak I",
    linewidth = 2,
    title = "Peak infected vs beta",
    label = "Peak I",
)

function plot_parameter_heatmap(df::DataFrame; metric::Symbol = :peak_I)
    beta_values = sort(unique(df.beta))
    gamma_values = sort(unique(df.gamma))
    matrix = fill(NaN, length(gamma_values), length(beta_values))

    for row in eachrow(df)
        i = findfirst(==(row.gamma), gamma_values)
        j = findfirst(==(row.beta), beta_values)
        matrix[i, j] = row[metric]
    end

    return heatmap(
        beta_values,
        gamma_values,
        matrix;
        xlabel = "beta",
        ylabel = "gamma",
        title = "Parameter grid: $(metric)",
        colorbar_title = string(metric),
    )
end

function sample_path(df::DataFrame, time_points)
    values = similar(collect(Float64.(time_points)))
    source_times = collect(df.time)
    for (idx, t) in enumerate(time_points)
        pos = searchsortedlast(source_times, t)
        pos = clamp(pos, 1, nrow(df))
        values[idx] = df.I[pos]
    end
    return values
end

function plot_comparison(df_det::DataFrame, df_stoch::DataFrame)
    sampled_I = sample_path(df_stoch, df_det.time)
    return plot(
        df_det.time,
        hcat(df_det.I, sampled_I);
        label = ["Deterministic I" "Stochastic I"],
        xlabel = "Time",
        ylabel = "Infected",
        linewidth = 2,
        title = "Deterministic vs stochastic trajectory",
    )
end

function plot_sir_network(net::LabelledPetriNet)
    xs = [-1.0, 0.0, 1.0]
    ys = [0.8, 0.8, 0.8]
    tx = [-0.5, 0.5]
    ty = [-0.4, -0.4]

    p = plot(
        xlim = (-1.4, 1.4),
        ylim = (-0.8, 1.2),
        legend = false,
        axis = nothing,
        grid = false,
        title = "Petri net for SIR",
        size = (700, 400),
    )

    scatter!(p, xs, ys; markersize = 14, marker = :circle, color = :steelblue)
    scatter!(p, tx, ty; markersize = 10, marker = :rect, color = :darkorange)

    annotations = [
        (-1.0, 0.8, text("S", 12, :white)),
        (0.0, 0.8, text("I", 12, :white)),
        (1.0, 0.8, text("R", 12, :white)),
        (-0.5, -0.4, text("infection", 10, :white)),
        (0.5, -0.4, text("recovery", 10, :white)),
    ]
    annotate!(p, annotations)

    plot!(p, [-1.0, -0.5], [0.68, -0.26]; color = :black, linewidth = 2)
    plot!(p, [0.0, -0.5], [0.68, -0.26]; color = :black, linewidth = 2)
    plot!(p, [-0.5, 0.0], [-0.26, 0.68]; color = :black, linewidth = 2)
    plot!(p, [-0.5, 0.08], [-0.26, 0.68]; color = :black, linewidth = 2)
    plot!(p, [0.0, 0.5], [0.68, -0.26]; color = :black, linewidth = 2)
    plot!(p, [0.5, 1.0], [-0.26, 0.68]; color = :black, linewidth = 2)

    annotate!(p, -0.72, 0.08, text("beta", 10))
    annotate!(p, 0.32, 0.08, text("gamma", 10))

    return p
end

function to_graphviz_sir(net::LabelledPetriNet)
    return """
digraph SIR {
  rankdir=LR;
  node [shape=circle, style=filled, fillcolor=lightblue];
  S;
  I;
  R;
  node [shape=box, style=filled, fillcolor=orange];
  infection [label="infection\\nbeta"];
  recovery [label="recovery\\ngamma"];
  S -> infection;
  I -> infection;
  infection -> I;
  infection -> I;
  I -> recovery;
  recovery -> R;
}
"""
end

function animate_sir(df::DataFrame, output_path::AbstractString; fps::Int = 8)
    anim = @animate for row in eachrow(df)
        bar(
            ["S", "I", "R"],
            [row.S, row.I, row.R];
            ylim = (0, maximum([maximum(df.S), maximum(df.I), maximum(df.R)]) * 1.05),
            color = [:steelblue, :firebrick, :seagreen],
            legend = false,
            ylabel = "Population",
            title = "t = $(round(row.time, digits = 1))",
        )
    end
    gif(anim, output_path; fps = fps)
    return output_path
end

end
