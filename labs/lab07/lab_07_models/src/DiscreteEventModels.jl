module DiscreteEventModels

using DataFrames
using Plots
using Random
using Statistics

export MMcParameters,
    RossParameters,
    compare_mmc_metrics,
    erlang_c_metrics,
    mmc_server_sweep,
    mmc_step_series,
    plot_healthy_machines,
    plot_mmc_comparison,
    plot_mmc_server_utilization,
    plot_mmc_sweep,
    plot_mmc_timeline,
    plot_mmc_wait_histogram,
    plot_ross_comparison,
    plot_ross_crash_distribution,
    plot_ross_queue,
    plot_ross_utilization_sweep,
    ross_machine_sweep,
    ross_mean_time_to_failure_analytic,
    ross_repairer_sweep,
    ross_replicates,
    simulate_mmc,
    simulate_ross,
    summarize_mmc

Base.@kwdef struct MMcParameters
    lambda::Float64 = 0.9
    mu::Float64 = 0.5
    servers::Int = 2
    num_customers::Int = 500
end

Base.@kwdef struct RossParameters
    machines::Int = 10
    spares::Int = 3
    repairers::Int = 1
    failure_mean::Float64 = 100.0
    repair_mean::Float64 = 1.0
end

sample_exponential(rng::AbstractRNG, mean_value::Real) = randexp(rng) * float(mean_value)

function validate_mmc(params::MMcParameters)
    params.lambda > 0 || throw(ArgumentError("lambda must be positive"))
    params.mu > 0 || throw(ArgumentError("mu must be positive"))
    params.servers > 0 || throw(ArgumentError("servers must be positive"))
    params.num_customers > 0 || throw(ArgumentError("num_customers must be positive"))
    params.lambda < params.servers * params.mu ||
        throw(ArgumentError("M/M/c requires lambda < c * mu for a stable regime"))
    return params
end

function validate_ross(params::RossParameters)
    params.machines > 0 || throw(ArgumentError("machines must be positive"))
    params.spares >= 0 || throw(ArgumentError("spares must be non-negative"))
    params.repairers > 0 || throw(ArgumentError("repairers must be positive"))
    params.failure_mean > 0 || throw(ArgumentError("failure_mean must be positive"))
    params.repair_mean > 0 || throw(ArgumentError("repair_mean must be positive"))
    return params
end

function erlang_c_metrics(lambda::Real, mu::Real, servers::Integer)
    params = validate_mmc(
        MMcParameters(
            lambda = float(lambda),
            mu = float(mu),
            servers = Int(servers),
            num_customers = 1,
        ),
    )

    ρ = params.lambda / (params.servers * params.mu)
    a = params.lambda / params.mu
    partial_sum = 1.0
    term = 1.0

    for n in 1:params.servers-1
        term *= a / n
        partial_sum += term
    end

    tail_base = term * a / params.servers
    tail = tail_base / (1.0 - ρ)
    p0 = 1.0 / (partial_sum + tail)
    p_wait = tail * p0
    Lq = p_wait * ρ / (1.0 - ρ)
    Wq = Lq / params.lambda
    W = Wq + 1.0 / params.mu
    L = params.lambda * W

    return (
        rho = ρ,
        p0 = p0,
        p_wait = p_wait,
        Lq = Lq,
        Wq = Wq,
        W = W,
        L = L,
    )
end

erlang_c_metrics(params::MMcParameters) = erlang_c_metrics(params.lambda, params.mu, params.servers)

function simulate_mmc(params::MMcParameters; rng::AbstractRNG = MersenneTwister(123))
    validate_mmc(params)

    server_available = zeros(Float64, params.servers)
    rows = NamedTuple[]
    arrival_time = 0.0

    for customer_id in 1:params.num_customers
        arrival_time += sample_exponential(rng, 1.0 / params.lambda)
        service_time = sample_exponential(rng, 1.0 / params.mu)
        server_id = argmin(server_available)
        service_start = max(arrival_time, server_available[server_id])
        departure_time = service_start + service_time
        waiting_time = service_start - arrival_time

        server_available[server_id] = departure_time

        push!(
            rows,
            (
                customer_id = customer_id,
                server_id = server_id,
                arrival_time = arrival_time,
                service_start = service_start,
                service_time = service_time,
                departure_time = departure_time,
                waiting_time = waiting_time,
                system_time = departure_time - arrival_time,
            ),
        )
    end

    return DataFrame(rows)
end

function events_to_step_dataframe(events; start_time::Real = 0.0)
    sort!(events, by = event -> (event.time, event.priority))
    times = Float64[float(start_time)]
    values = Int[0]
    current = 0
    idx = 1

    while idx <= length(events)
        event_time = events[idx].time
        while idx <= length(events) && abs(events[idx].time - event_time) <= 1.0e-12
            current += events[idx].delta
            idx += 1
        end
        push!(times, event_time)
        push!(values, current)
    end

    return DataFrame(time = times, value = values)
end

function mmc_step_series(df::DataFrame)
    queue_events = NamedTuple[]
    service_events = NamedTuple[]
    system_events = NamedTuple[]

    for row in eachrow(df)
        push!(system_events, (time = row.arrival_time, delta = 1, priority = 2))
        push!(system_events, (time = row.departure_time, delta = -1, priority = 1))

        push!(service_events, (time = row.service_start, delta = 1, priority = 2))
        push!(service_events, (time = row.departure_time, delta = -1, priority = 1))

        if row.waiting_time > 1.0e-12
            push!(queue_events, (time = row.arrival_time, delta = 1, priority = 2))
            push!(queue_events, (time = row.service_start, delta = -1, priority = 1))
        end
    end

    start_time = minimum(df.arrival_time)
    return (
        queue = events_to_step_dataframe(queue_events; start_time = start_time),
        service = events_to_step_dataframe(service_events; start_time = start_time),
        system = events_to_step_dataframe(system_events; start_time = start_time),
    )
end

function summarize_mmc(df::DataFrame, params::MMcParameters)
    validate_mmc(params)

    horizon = maximum(df.departure_time)
    avg_waiting_time = mean(df.waiting_time)
    avg_system_time = mean(df.system_time)
    wait_probability = mean(df.waiting_time .> 1.0e-12)
    avg_queue_length = params.lambda * avg_waiting_time
    avg_system_size = params.lambda * avg_system_time

    busy_time_by_server = zeros(Float64, params.servers)
    for row in eachrow(df)
        busy_time_by_server[row.server_id] += row.service_time
    end

    server_utilization = busy_time_by_server ./ horizon

    return (
        avg_waiting_time = avg_waiting_time,
        avg_system_time = avg_system_time,
        wait_probability = wait_probability,
        avg_queue_length = avg_queue_length,
        avg_system_size = avg_system_size,
        horizon = horizon,
        mean_server_utilization = mean(server_utilization),
    )
end

function compare_mmc_metrics(df::DataFrame, params::MMcParameters)
    sim = summarize_mmc(df, params)
    analytic = erlang_c_metrics(params)

    return DataFrame(
        metric = ["P(wait)", "Wq", "W", "Lq", "L"],
        simulated = [
            sim.wait_probability,
            sim.avg_waiting_time,
            sim.avg_system_time,
            sim.avg_queue_length,
            sim.avg_system_size,
        ],
        analytic = [
            analytic.p_wait,
            analytic.Wq,
            analytic.W,
            analytic.Lq,
            analytic.L,
        ],
    )
end

function plot_mmc_timeline(df::DataFrame; title::AbstractString = "M/M/c timeline")
    series = mmc_step_series(df)
    plt = plot(
        series.system.time,
        series.system.value;
        seriestype = :steppost,
        label = "In system",
        linewidth = 2,
        xlabel = "Time",
        ylabel = "Customers",
        title = title,
    )
    plot!(
        plt,
        series.queue.time,
        series.queue.value;
        seriestype = :steppost,
        label = "In queue",
        linewidth = 2,
    )
    plot!(
        plt,
        series.service.time,
        series.service.value;
        seriestype = :steppost,
        label = "In service",
        linewidth = 2,
    )
    return plt
end

function plot_mmc_wait_histogram(df::DataFrame; title::AbstractString = "Waiting time distribution")
    return histogram(
        df.waiting_time;
        bins = :sturges,
        xlabel = "Waiting time",
        ylabel = "Frequency",
        title = title,
        label = "Customers",
        normalize = false,
    )
end

function plot_mmc_server_utilization(
    df::DataFrame,
    params::MMcParameters;
    title::AbstractString = "Server utilization",
)
    summary = summarize_mmc(df, params)
    busy_time_by_server = zeros(Float64, params.servers)
    for row in eachrow(df)
        busy_time_by_server[row.server_id] += row.service_time
    end
    utilization = busy_time_by_server ./ summary.horizon

    return bar(
        1:params.servers,
        utilization;
        xlabel = "Server",
        ylabel = "Utilization",
        title = title,
        label = "Busy fraction",
        ylim = (0.0, 1.05),
    )
end

function plot_mmc_comparison(df::DataFrame, params::MMcParameters)
    comparison = compare_mmc_metrics(df, params)
    return bar(
        comparison.metric,
        hcat(comparison.simulated, comparison.analytic);
        label = ["Simulation" "Analytic"],
        xlabel = "Metric",
        ylabel = "Value",
        title = "M/M/c: simulation vs analytic formula",
    )
end

function mmc_server_sweep(
    server_counts;
    lambda::Real = 0.9,
    mu::Real = 0.5,
    num_customers::Int = 800,
    seed::Integer = 123,
)
    rows = NamedTuple[]

    for (idx, servers) in enumerate(server_counts)
        params = MMcParameters(
            lambda = float(lambda),
            mu = float(mu),
            servers = Int(servers),
            num_customers = num_customers,
        )
        df = simulate_mmc(params; rng = MersenneTwister(seed + idx - 1))
        sim = summarize_mmc(df, params)
        analytic = erlang_c_metrics(params)

        push!(
            rows,
            (
                servers = params.servers,
                sim_wait_probability = sim.wait_probability,
                analytic_wait_probability = analytic.p_wait,
                sim_avg_waiting = sim.avg_waiting_time,
                analytic_avg_waiting = analytic.Wq,
                sim_avg_system = sim.avg_system_time,
                analytic_avg_system = analytic.W,
            ),
        )
    end

    return DataFrame(rows)
end

function plot_mmc_sweep(df::DataFrame)
    p1 = plot(
        df.servers,
        hcat(df.sim_wait_probability, df.analytic_wait_probability);
        label = ["Simulation" "Analytic"],
        linewidth = 2,
        marker = :circle,
        xlabel = "Servers",
        ylabel = "P(wait)",
        title = "Waiting probability vs number of servers",
    )
    p2 = plot(
        df.servers,
        hcat(df.sim_avg_waiting, df.analytic_avg_waiting);
        label = ["Simulation" "Analytic"],
        linewidth = 2,
        marker = :square,
        xlabel = "Servers",
        ylabel = "Average waiting time",
        title = "Average queueing delay vs number of servers",
    )
    return plot(p1, p2; layout = (2, 1), size = (850, 700))
end

function ross_mean_time_to_failure_analytic(params::RossParameters)
    validate_ross(params)

    n_states = params.spares + 1
    A = zeros(Float64, n_states, n_states)
    b = ones(Float64, n_states)
    top_state = params.machines + params.spares

    for row in 1:n_states
        healthy = params.machines + row - 1
        failures = min(params.machines, healthy) / params.failure_mean
        broken = top_state - healthy
        repairs = min(params.repairers, broken) / params.repair_mean

        A[row, row] = failures + repairs
        row > 1 && (A[row, row - 1] = -failures)
        row < n_states && (A[row, row + 1] = -repairs)
    end

    expected_times = A \ b
    return expected_times[end]
end

function simulate_ross(params::RossParameters; seed::Integer = 123)
    validate_ross(params)

    rng = MersenneTwister(seed)
    time = 0.0
    healthy = params.machines + params.spares
    spares = params.spares
    busy_repairers = 0
    repair_queue = 0
    repair_completion_times = Float64[]

    area_busy = 0.0
    area_queue = 0.0
    area_healthy = 0.0

    rows = NamedTuple[
        (
            time = 0.0,
            healthy_machines = healthy,
            working_machines = params.machines,
            spare_machines = spares,
            busy_repairers = busy_repairers,
            repair_queue = repair_queue,
            event = "start",
        ),
    ]

    while true
        failure_rate = params.machines / params.failure_mean
        next_failure_time = time + sample_exponential(rng, 1.0 / failure_rate)
        next_repair_time = isempty(repair_completion_times) ? Inf : minimum(repair_completion_times)
        next_time = min(next_failure_time, next_repair_time)
        dt = next_time - time

        area_busy += busy_repairers * dt
        area_queue += repair_queue * dt
        area_healthy += healthy * dt

        time = next_time

        if next_failure_time <= next_repair_time
            if spares == 0
                healthy -= 1
                push!(
                    rows,
                    (
                        time = time,
                        healthy_machines = healthy,
                        working_machines = min(params.machines, healthy),
                        spare_machines = max(healthy - params.machines, 0),
                        busy_repairers = busy_repairers,
                        repair_queue = repair_queue,
                        event = "crash",
                    ),
                )
                break
            end

            spares -= 1
            healthy -= 1
            if busy_repairers < params.repairers
                busy_repairers += 1
                push!(repair_completion_times, time + sample_exponential(rng, params.repair_mean))
            else
                repair_queue += 1
            end

            push!(
                rows,
                (
                    time = time,
                    healthy_machines = healthy,
                    working_machines = params.machines,
                    spare_machines = spares,
                    busy_repairers = busy_repairers,
                    repair_queue = repair_queue,
                    event = "failure",
                ),
            )
        else
            completion_index = argmin(repair_completion_times)
            deleteat!(repair_completion_times, completion_index)

            healthy += 1
            spares += 1
            if repair_queue > 0
                repair_queue -= 1
                push!(repair_completion_times, time + sample_exponential(rng, params.repair_mean))
            else
                busy_repairers -= 1
            end

            push!(
                rows,
                (
                    time = time,
                    healthy_machines = healthy,
                    working_machines = params.machines,
                    spare_machines = spares,
                    busy_repairers = busy_repairers,
                    repair_queue = repair_queue,
                    event = "repair_complete",
                ),
            )
        end
    end

    return (
        trace = DataFrame(rows),
        crash_time = time,
        utilization = area_busy / (params.repairers * time),
        mean_queue_length = area_queue / time,
        mean_healthy = area_healthy / time,
    )
end

function ross_replicates(
    params::RossParameters;
    runs::Int = 50,
    seed::Integer = 123,
)
    runs > 0 || throw(ArgumentError("runs must be positive"))
    rows = NamedTuple[]

    for run in 1:runs
        result = simulate_ross(params; seed = seed + run - 1)
        push!(
            rows,
            (
                run = run,
                crash_time = result.crash_time,
                utilization = result.utilization,
                mean_queue_length = result.mean_queue_length,
                mean_healthy = result.mean_healthy,
            ),
        )
    end

    return DataFrame(rows)
end

function ross_machine_sweep(
    machine_counts;
    spares::Int = 3,
    repairers::Int = 2,
    failure_mean::Real = 100.0,
    repair_mean::Real = 1.0,
    runs::Int = 40,
    seed::Integer = 123,
)
    rows = NamedTuple[]

    for (idx, machines) in enumerate(machine_counts)
        params = RossParameters(
            machines = Int(machines),
            spares = spares,
            repairers = repairers,
            failure_mean = float(failure_mean),
            repair_mean = float(repair_mean),
        )
        replications = ross_replicates(params; runs = runs, seed = seed + 100 * (idx - 1))
        push!(
            rows,
            (
                machines = params.machines,
                sim_mean_crash_time = mean(replications.crash_time),
                sim_std_crash_time = std(replications.crash_time),
                analytic_crash_time = ross_mean_time_to_failure_analytic(params),
                mean_utilization = mean(replications.utilization),
                mean_queue_length = mean(replications.mean_queue_length),
            ),
        )
    end

    return DataFrame(rows)
end

function ross_repairer_sweep(
    repairer_counts;
    machines::Int = 10,
    spares::Int = 3,
    failure_mean::Real = 100.0,
    repair_mean::Real = 1.0,
    runs::Int = 40,
    seed::Integer = 123,
)
    rows = NamedTuple[]

    for (idx, repairers) in enumerate(repairer_counts)
        params = RossParameters(
            machines = machines,
            spares = spares,
            repairers = Int(repairers),
            failure_mean = float(failure_mean),
            repair_mean = float(repair_mean),
        )
        replications = ross_replicates(params; runs = runs, seed = seed + 100 * (idx - 1))
        push!(
            rows,
            (
                repairers = params.repairers,
                sim_mean_crash_time = mean(replications.crash_time),
                sim_std_crash_time = std(replications.crash_time),
                analytic_crash_time = ross_mean_time_to_failure_analytic(params),
                mean_utilization = mean(replications.utilization),
                mean_queue_length = mean(replications.mean_queue_length),
            ),
        )
    end

    return DataFrame(rows)
end

function plot_healthy_machines(trace::DataFrame; title::AbstractString = "Healthy machines over time")
    return plot(
        trace.time,
        hcat(trace.healthy_machines, trace.working_machines);
        seriestype = :steppost,
        label = ["Healthy" "Working"],
        linewidth = 2,
        xlabel = "Time",
        ylabel = "Machines",
        title = title,
    )
end

function plot_ross_queue(trace::DataFrame; title::AbstractString = "Repair facility load")
    return plot(
        trace.time,
        hcat(trace.busy_repairers, trace.repair_queue);
        seriestype = :steppost,
        label = ["Busy repairers" "Queue length"],
        linewidth = 2,
        xlabel = "Time",
        ylabel = "Count",
        title = title,
    )
end

function plot_ross_crash_distribution(
    replications::DataFrame;
    title::AbstractString = "Crash-time distribution",
)
    return histogram(
        replications.crash_time;
        bins = :sturges,
        xlabel = "Crash time",
        ylabel = "Frequency",
        title = title,
        label = "Runs",
    )
end

function plot_ross_comparison(
    df::DataFrame;
    xcol::Symbol = :machines,
    title::AbstractString = "Ross model: mean crash time",
)
    x = df[!, xcol]
    return plot(
        x,
        hcat(df.sim_mean_crash_time, df.analytic_crash_time);
        linewidth = 2,
        marker = [:circle :diamond],
        label = ["Simulation" "Analytic"],
        xlabel = string(xcol),
        ylabel = "Mean time to crash",
        title = title,
    )
end

function plot_ross_utilization_sweep(
    df::DataFrame;
    xcol::Symbol = :machines,
    title::AbstractString = "Ross model: repair load",
)
    p1 = plot(
        df[!, xcol],
        df.mean_utilization;
        linewidth = 2,
        marker = :circle,
        label = "Utilization",
        xlabel = string(xcol),
        ylabel = "Busy fraction",
        title = "Repairer utilization",
    )
    p2 = plot(
        df[!, xcol],
        df.mean_queue_length;
        linewidth = 2,
        marker = :square,
        label = "Queue length",
        xlabel = string(xcol),
        ylabel = "Average queue",
        title = "Average repair queue length",
    )
    return plot(p1, p2; layout = (2, 1), size = (850, 700), plot_title = title)
end

end
