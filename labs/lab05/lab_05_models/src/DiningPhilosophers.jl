module DiningPhilosophers

using DataFrames
using Plots
using Random
using Statistics

export PetriNet
export build_classical_network, build_arbiter_network
export enabled_transitions, simulate_stochastic, simulate_deterministic
export detect_deadlock, deadlock_time, total_eaters
export plot_marking_evolution, plot_network, plot_eating_comparison
export plot_total_eaters_comparison, animate_marking
export parameter_scan, summarize_parameter_scan, plot_parameter_summary

struct PetriNet
    n_places::Int
    n_transitions::Int
    pre::Matrix{Int}
    post::Matrix{Int}
    incidence::Matrix{Int}
    place_names::Vector{Symbol}
    transition_names::Vector{Symbol}
    metadata::Dict{Symbol, Any}
end

function PetriNet(
    n_places::Integer,
    n_transitions::Integer;
    place_names = Symbol[],
    transition_names = Symbol[],
    metadata = Dict{Symbol, Any}(),
)
    names_places = isempty(place_names) ? [Symbol("p$i") for i in 1:n_places] : copy(place_names)
    names_transitions = isempty(transition_names) ? [Symbol("t$i") for i in 1:n_transitions] : copy(transition_names)
    pre = zeros(Int, n_places, n_transitions)
    post = zeros(Int, n_places, n_transitions)
    incidence = zeros(Int, n_places, n_transitions)
    return PetriNet(
        Int(n_places),
        Int(n_transitions),
        pre,
        post,
        incidence,
        names_places,
        names_transitions,
        Dict{Symbol, Any}(metadata),
    )
end

function add_input_arc!(net::PetriNet, place::Int, transition::Int, weight::Int = 1)
    net.pre[place, transition] += weight
    net.incidence[place, transition] -= weight
    return net
end

function add_output_arc!(net::PetriNet, place::Int, transition::Int, weight::Int = 1)
    net.post[place, transition] += weight
    net.incidence[place, transition] += weight
    return net
end

function build_classical_network(N::Int)
    n_places = 4 * N
    n_transitions = 3 * N
    net = PetriNet(
        n_places,
        n_transitions;
        metadata = Dict(:N => N, :model => :classic, :has_arbiter => false),
    )

    for i in 1:N
        net.place_names[i] = Symbol("Think_$i")
        net.place_names[N + i] = Symbol("Hungry_$i")
        net.place_names[2N + i] = Symbol("Eat_$i")
        net.place_names[3N + i] = Symbol("Fork_$i")

        net.transition_names[i] = Symbol("GetLeft_$i")
        net.transition_names[N + i] = Symbol("GetRight_$i")
        net.transition_names[2N + i] = Symbol("PutForks_$i")
    end

    for i in 1:N
        think = i
        hungry = N + i
        eat = 2N + i
        left_fork = 3N + i
        right_fork = 3N + (i % N + 1)
        get_left = i
        get_right = N + i
        put_forks = 2N + i

        add_input_arc!(net, think, get_left)
        add_input_arc!(net, left_fork, get_left)
        add_output_arc!(net, hungry, get_left)

        add_input_arc!(net, hungry, get_right)
        add_input_arc!(net, right_fork, get_right)
        add_output_arc!(net, eat, get_right)

        add_input_arc!(net, eat, put_forks)
        add_output_arc!(net, think, put_forks)
        add_output_arc!(net, left_fork, put_forks)
        add_output_arc!(net, right_fork, put_forks)
    end

    u0 = zeros(Float64, n_places)
    for i in 1:N
        u0[i] = 1.0
        u0[3N + i] = 1.0
    end

    return net, u0, copy(net.place_names)
end

function build_arbiter_network(N::Int)
    n_places = 4 * N + 1
    n_transitions = 3 * N
    net = PetriNet(
        n_places,
        n_transitions;
        metadata = Dict(:N => N, :model => :arbiter, :has_arbiter => true),
    )

    for i in 1:N
        net.place_names[i] = Symbol("Think_$i")
        net.place_names[N + i] = Symbol("Hungry_$i")
        net.place_names[2N + i] = Symbol("Eat_$i")
        net.place_names[3N + i] = Symbol("Fork_$i")

        net.transition_names[i] = Symbol("GetLeft_$i")
        net.transition_names[N + i] = Symbol("GetRight_$i")
        net.transition_names[2N + i] = Symbol("PutForks_$i")
    end
    net.place_names[4N + 1] = :Arbiter

    arbiter_idx = 4N + 1
    for i in 1:N
        think = i
        hungry = N + i
        eat = 2N + i
        left_fork = 3N + i
        right_fork = 3N + (i % N + 1)
        get_left = i
        get_right = N + i
        put_forks = 2N + i

        add_input_arc!(net, think, get_left)
        add_input_arc!(net, left_fork, get_left)
        add_input_arc!(net, arbiter_idx, get_left)
        add_output_arc!(net, hungry, get_left)

        add_input_arc!(net, hungry, get_right)
        add_input_arc!(net, right_fork, get_right)
        add_output_arc!(net, eat, get_right)

        add_input_arc!(net, eat, put_forks)
        add_output_arc!(net, think, put_forks)
        add_output_arc!(net, left_fork, put_forks)
        add_output_arc!(net, right_fork, put_forks)
        add_output_arc!(net, arbiter_idx, put_forks)
    end

    u0 = zeros(Float64, n_places)
    for i in 1:N
        u0[i] = 1.0
        u0[3N + i] = 1.0
    end
    u0[arbiter_idx] = N - 1

    return net, u0, copy(net.place_names)
end

function enabled_transitions(net::PetriNet, marking::AbstractVector{<:Real}; tol::Real = 1.0e-9)
    enabled = Int[]
    for j in 1:net.n_transitions
        can_fire = true
        for i in 1:net.n_places
            if net.pre[i, j] > 0 && marking[i] + tol < net.pre[i, j]
                can_fire = false
                break
            end
        end
        can_fire && push!(enabled, j)
    end
    return enabled
end

function transition_propensities(
    net::PetriNet,
    marking::AbstractVector{<:Real},
    rates::AbstractVector{<:Real},
)
    a = zeros(Float64, net.n_transitions)
    for j in 1:net.n_transitions
        prod = float(rates[j])
        for i in 1:net.n_places
            if net.pre[i, j] > 0
                available = max(marking[i], 0.0)
                if available + 1.0e-9 < net.pre[i, j]
                    prod = 0.0
                    break
                end
                prod *= available ^ net.pre[i, j]
            end
        end
        a[j] = prod
    end
    return a
end

function fire_transition!(net::PetriNet, marking::Vector{Float64}, transition::Int)
    marking .+= @view net.incidence[:, transition]
    return marking
end

function marking_dataframe(times, states, place_names)
    df = DataFrame(time = times)
    for (idx, name) in enumerate(place_names)
        df[!, String(name)] = [state[idx] for state in states]
    end
    return df
end

function simulate_stochastic(
    net::PetriNet,
    u0::AbstractVector{<:Real},
    tmax::Real;
    rates = ones(Float64, net.n_transitions),
    rng = Random.default_rng(),
    max_steps::Int = 200_000,
)
    marking = Float64.(u0)
    t = 0.0
    times = [t]
    states = [copy(marking)]
    steps = 0

    while t < tmax && steps < max_steps
        steps += 1
        a = transition_propensities(net, marking, rates)
        a0 = sum(a)
        if a0 <= 1.0e-12
            break
        end

        dt = -log(rand(rng)) / a0
        threshold = rand(rng) * a0
        cumsum = 0.0
        chosen = net.n_transitions
        for j in 1:net.n_transitions
            cumsum += a[j]
            if threshold <= cumsum
                chosen = j
                break
            end
        end

        fire_transition!(net, marking, chosen)
        t += dt
        if t <= tmax + 1.0e-9
            push!(times, min(t, float(tmax)))
            push!(states, copy(marking))
        end
    end

    return marking_dataframe(times, states, net.place_names)
end

function simulate_deterministic(
    net::PetriNet,
    u0::AbstractVector{<:Real},
    tmax::Real;
    rates = ones(Float64, net.n_transitions),
    saveat::Real = 0.2,
    dt::Real = 0.02,
)
    save_times = collect(0.0:saveat:float(tmax))
    if isempty(save_times) || save_times[end] < tmax
        push!(save_times, float(tmax))
    end

    marking = Float64.(u0)
    times = Float64[0.0]
    states = [copy(marking)]
    t = 0.0

    for target in save_times[2:end]
        while t + 1.0e-12 < target
            h = min(float(dt), target - t)
            a = transition_propensities(net, marking, rates)
            derivative = net.incidence * a
            marking .= max.(marking .+ h .* derivative, 0.0)
            t += h
        end
        push!(times, target)
        push!(states, copy(marking))
    end

    return marking_dataframe(times, states, net.place_names)
end

function last_marking(df::DataFrame, net::PetriNet)
    return [df[end, String(net.place_names[i])] for i in 1:net.n_places]
end

function detect_deadlock(net::PetriNet, marking::AbstractVector{<:Real}; tol::Real = 1.0e-9)
    return isempty(enabled_transitions(net, marking; tol = tol))
end

function detect_deadlock(df::DataFrame, net::PetriNet; tol::Real = 1.0e-9)
    return detect_deadlock(net, last_marking(df, net); tol = tol)
end

function deadlock_time(df::DataFrame, net::PetriNet)
    return detect_deadlock(df, net) ? df[end, :time] : missing
end

function total_eaters(df::DataFrame, N::Int)
    cols = [Symbol("Eat_$i") for i in 1:N]
    return vec(sum(Matrix(df[:, cols]); dims = 2))
end

function plot_marking_evolution(df::DataFrame, net::PetriNet)
    N = Int(net.metadata[:N])
    groups = [
        ("Think", 1:N),
        ("Hungry", N + 1:2N),
        ("Eat", 2N + 1:3N),
        ("Fork", 3N + 1:4N),
    ]
    get(net.metadata, :has_arbiter, false) && push!(groups, ("Arbiter", 4N + 1:4N + 1))

    panels = Plots.Plot[]
    for (group, range) in groups
        p = plot(
            xlabel = "Время",
            ylabel = "Фишки",
            title = group,
            legend = :right,
        )
        for idx in range
            col = String(net.place_names[idx])
            plot!(p, df.time, df[!, col], label = col, lw = 2)
        end
        push!(panels, p)
    end

    return plot(panels...; layout = (length(groups), 1), size = (900, 240 * length(groups)))
end

function plot_eating_comparison(df_classic::DataFrame, df_arbiter::DataFrame, N::Int)
    eat_cols = [Symbol("Eat_$i") for i in 1:N]
    labels = ["Ф$i" for i in 1:N]

    p1 = plot(
        df_classic.time,
        Matrix(df_classic[:, eat_cols]);
        label = labels,
        xlabel = "Время",
        ylabel = "Eat_i",
        title = "Классическая сеть",
        lw = 2,
    )
    p2 = plot(
        df_arbiter.time,
        Matrix(df_arbiter[:, eat_cols]);
        label = labels,
        xlabel = "Время",
        ylabel = "Eat_i",
        title = "Сеть с арбитром",
        lw = 2,
    )

    return plot(p1, p2; layout = (2, 1), size = (900, 700))
end

function plot_total_eaters_comparison(
    df_classic_stochastic::DataFrame,
    df_classic_deterministic::DataFrame,
    df_arbiter_stochastic::DataFrame,
    df_arbiter_deterministic::DataFrame,
    N::Int,
)
    p1 = plot(
        df_classic_stochastic.time,
        total_eaters(df_classic_stochastic, N);
        label = "Стохастическая",
        xlabel = "Время",
        ylabel = "Сумма Eat_i",
        title = "Классическая сеть",
        lw = 2,
    )
    plot!(p1, df_classic_deterministic.time, total_eaters(df_classic_deterministic, N); label = "Детерминированная", lw = 2, ls = :dash)

    p2 = plot(
        df_arbiter_stochastic.time,
        total_eaters(df_arbiter_stochastic, N);
        label = "Стохастическая",
        xlabel = "Время",
        ylabel = "Сумма Eat_i",
        title = "Сеть с арбитром",
        lw = 2,
    )
    plot!(p2, df_arbiter_deterministic.time, total_eaters(df_arbiter_deterministic, N); label = "Детерминированная", lw = 2, ls = :dash)

    return plot(p1, p2; layout = (2, 1), size = (900, 650))
end

function plot_network(net::PetriNet)
    N = Int(net.metadata[:N])
    has_arbiter = get(net.metadata, :has_arbiter, false)

    place_x = Float64[]
    place_y = Float64[]
    place_label = String[]
    place_color = Symbol[]

    transition_x = Float64[]
    transition_y = Float64[]
    transition_label = String[]

    coords = Dict{Symbol, Tuple{Float64, Float64}}()

    for i in 1:N
        coords[Symbol("Think_$i")] = (i, 3.3)
        coords[Symbol("Hungry_$i")] = (i, 2.2)
        coords[Symbol("Eat_$i")] = (i, 1.1)
        coords[Symbol("Fork_$i")] = (i, 0.0)

        coords[Symbol("GetLeft_$i")] = (i - 0.22, 2.75)
        coords[Symbol("GetRight_$i")] = (i + 0.22, 1.65)
        coords[Symbol("PutForks_$i")] = (i, 0.55)
    end
    if has_arbiter
        coords[:Arbiter] = ((N + 1) / 2, 4.25)
    end

    edges = Tuple{Symbol, Symbol}[]
    for i in 1:N
        right_fork = Symbol("Fork_$(i % N + 1)")
        push!(edges, (Symbol("Think_$i"), Symbol("GetLeft_$i")))
        push!(edges, (Symbol("Fork_$i"), Symbol("GetLeft_$i")))
        push!(edges, (Symbol("GetLeft_$i"), Symbol("Hungry_$i")))
        push!(edges, (Symbol("Hungry_$i"), Symbol("GetRight_$i")))
        push!(edges, (right_fork, Symbol("GetRight_$i")))
        push!(edges, (Symbol("GetRight_$i"), Symbol("Eat_$i")))
        push!(edges, (Symbol("Eat_$i"), Symbol("PutForks_$i")))
        push!(edges, (Symbol("PutForks_$i"), Symbol("Think_$i")))
        push!(edges, (Symbol("PutForks_$i"), Symbol("Fork_$i")))
        push!(edges, (Symbol("PutForks_$i"), right_fork))

        if has_arbiter
            push!(edges, (:Arbiter, Symbol("GetLeft_$i")))
            push!(edges, (Symbol("PutForks_$i"), :Arbiter))
        end
    end

    p = plot(
        legend = false,
        xlim = (0.4, N + 0.6),
        ylim = (-0.5, has_arbiter ? 4.8 : 3.8),
        axis = false,
        title = has_arbiter ? "Сеть Петри с арбитром" : "Классическая сеть Петри",
        size = (1000, has_arbiter ? 700 : 620),
    )

    for (src, dst) in edges
        x1, y1 = coords[src]
        x2, y2 = coords[dst]
        plot!(p, [x1, x2], [y1, y2]; c = :gray55, lw = 1.7)
    end

    for i in 1:N
        for (name, color) in (
            (Symbol("Think_$i"), :steelblue),
            (Symbol("Hungry_$i"), :goldenrod1),
            (Symbol("Eat_$i"), :seagreen3),
            (Symbol("Fork_$i"), :indianred2),
        )
            x, y = coords[name]
            push!(place_x, x)
            push!(place_y, y)
            push!(place_label, String(name))
            push!(place_color, color)
        end

        for name in (Symbol("GetLeft_$i"), Symbol("GetRight_$i"), Symbol("PutForks_$i"))
            x, y = coords[name]
            push!(transition_x, x)
            push!(transition_y, y)
            push!(transition_label, String(name))
        end
    end

    scatter!(
        p,
        place_x,
        place_y;
        marker = (:circle, 11),
        markercolor = place_color,
        markerstrokecolor = :black,
    )
    scatter!(
        p,
        transition_x,
        transition_y;
        marker = (:rect, 10),
        markercolor = :white,
        markerstrokecolor = :black,
    )

    for (x, y, label) in zip(place_x, place_y, place_label)
        annotate!(p, x, y + 0.18, text(label, 8, :black))
    end
    for (x, y, label) in zip(transition_x, transition_y, transition_label)
        annotate!(p, x, y + 0.18, text(label, 7, :black))
    end

    if has_arbiter
        x, y = coords[:Arbiter]
        scatter!(p, [x], [y]; marker = (:circle, 12), markercolor = :mediumpurple1, markerstrokecolor = :black)
        annotate!(p, x, y + 0.2, text("Arbiter", 9, :black))
    end

    return p
end

function animate_marking(df::DataFrame, net::PetriNet, output_path::AbstractString; fps::Int = 2)
    labels = [String(name) for name in net.place_names]
    ymax = maximum(Matrix(df[:, Not(:time)])) + 1
    anim = @animate for row in eachrow(df)
        values = [row[Symbol(label)] for label in labels]
        bar(
            1:length(values),
            values;
            legend = false,
            xlabel = "Позиция",
            ylabel = "Фишки",
            ylims = (0, ymax),
            title = "Время = $(round(row.time; digits = 2))",
            color = :steelblue,
        )
        xticks!(1:length(values), labels, xrotation = 45)
    end
    gif(anim, output_path; fps = fps)
    return output_path
end

function parameter_scan(
    philosopher_counts;
    seeds = 1:20,
    tmax::Real = 60.0,
)
    rows = NamedTuple[]
    for N in philosopher_counts
        for (model_name, builder) in (("classic", build_classical_network), ("arbiter", build_arbiter_network))
            for seed in seeds
                net, u0, _ = builder(N)
                df = simulate_stochastic(net, u0, tmax; rng = MersenneTwister(seed))
                push!(
                    rows,
                    (
                        model = model_name,
                        philosophers = N,
                        seed = seed,
                        deadlock = detect_deadlock(df, net),
                        deadlock_time = deadlock_time(df, net),
                        final_eaters = total_eaters(df, N)[end],
                        final_hungry = sum(df[end, Symbol("Hungry_$i")] for i in 1:N),
                        events = nrow(df) - 1,
                    ),
                )
            end
        end
    end
    return DataFrame(rows)
end

function summarize_parameter_scan(raw::DataFrame)
    grouped = groupby(raw, [:model, :philosophers])
    return combine(
        grouped,
        :deadlock => mean => :deadlock_probability,
        :deadlock_time => (x -> isempty(skipmissing(x)) ? missing : mean(skipmissing(x))) => :mean_deadlock_time,
        :final_eaters => mean => :mean_final_eaters,
        :final_hungry => mean => :mean_final_hungry,
        :events => mean => :mean_events,
    )
end

function plot_parameter_summary(summary::DataFrame)
    classic = filter(row -> row.model == "classic", summary)
    arbiter = filter(row -> row.model == "arbiter", summary)

    p1 = plot(
        classic.philosophers,
        classic.deadlock_probability;
        label = "Классическая сеть",
        marker = :circle,
        lw = 2,
        xlabel = "Число философов N",
        ylabel = "Вероятность deadlock",
        title = "Вероятность взаимной блокировки",
        ylim = (-0.05, 1.05),
    )
    plot!(p1, arbiter.philosophers, arbiter.deadlock_probability; label = "С арбитром", marker = :square, lw = 2)

    classic_times = coalesce.(classic.mean_deadlock_time, NaN)
    arbiter_times = coalesce.(arbiter.mean_deadlock_time, NaN)
    p2 = plot(
        classic.philosophers,
        classic_times;
        label = "Классическая сеть",
        marker = :circle,
        lw = 2,
        xlabel = "Число философов N",
        ylabel = "Среднее время deadlock",
        title = "Среднее время наступления deadlock",
    )
    plot!(p2, arbiter.philosophers, arbiter_times; label = "С арбитром", marker = :square, lw = 2)

    return plot(p1, p2; layout = (2, 1), size = (900, 700))
end

end
