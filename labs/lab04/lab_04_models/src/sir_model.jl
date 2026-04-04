using Agents
using Agents.Graphs
using Distributions: Poisson
using DrWatson: @dict
using Random
using StatsBase: Weights, sample

@agent struct Person(GraphAgent)
    days_infected::Int
    status::Symbol
    home_city::Int
end

function default_migration_matrix(Ns)
    C = length(Ns)
    migration_rates = zeros(Float64, C, C)
    for i in 1:C
        for j in 1:C
            migration_rates[i, j] = (Ns[i] + Ns[j]) / Ns[i]
        end
        migration_rates[i, :] ./= sum(migration_rates[i, :])
    end
    return migration_rates
end

function create_migration_matrix(C::Integer, intensity::Real)
    @assert C >= 2 "At least two cities are required"
    @assert 0.0 <= intensity <= 1.0 "Migration intensity must be in [0, 1]"

    migration_rates = fill(intensity / (C - 1), C, C)
    for i in 1:C
        migration_rates[i, i] = 1 - intensity
    end
    return migration_rates
end

function initialize_sir(;
    Ns = [1000, 1000, 1000],
    migration_rates = nothing,
    beta_und = [0.5, 0.5, 0.5],
    beta_det = [0.05, 0.05, 0.05],
    infection_period = 14,
    detection_time = 7,
    death_rate = 0.02,
    reinfection_probability = 0.1,
    Is = [0, 0, 1],
    seed = 42,
    n_steps = 100,
    quarantine_threshold = nothing,
)
    C = length(Ns)
    @assert C == length(beta_und) == length(beta_det) == length(Is) "Parameter lengths must match"

    rng = Xoshiro(seed)
    migration_rates === nothing && (migration_rates = default_migration_matrix(Ns))
    migration_rates = Float64.(migration_rates)
    city_locked_down = falses(C)

    properties = @dict(
        Ns,
        beta_und,
        beta_det,
        migration_rates,
        infection_period,
        detection_time,
        death_rate,
        reinfection_probability,
        C,
        n_steps,
        quarantine_threshold,
        city_locked_down,
    )

    space = GraphSpace(complete_graph(C))
    model = StandardABM(Person, space; properties, rng, agent_step! = sir_agent_step!)

    for city in 1:C
        for _ in 1:Ns[city]
            add_agent!(city, model, 0, :S, city)
        end
    end

    for city in 1:C
        Is[city] == 0 && continue
        city_agents = ids_in_position(city, model)
        infected_ids = sample(rng, city_agents, Is[city]; replace = false)
        for id in infected_ids
            agent = model[id]
            agent.status = :I
            agent.days_infected = 1
        end
    end

    return model
end

function sir_agent_step!(agent, model)
    migrate!(agent, model)
    agent.status == :I && transmit!(agent, model)
    agent.status == :I && (agent.days_infected += 1)
    recover_or_die!(agent, model)
    return nothing
end

function migrate!(agent, model)
    current_city = agent.pos
    probs = model.migration_rates[current_city, :]
    target_city = sample(abmrng(model), 1:model.C, Weights(probs))
    target_city == current_city && return nothing
    move_agent!(agent, target_city, model)
    return nothing
end

function transmit!(agent, model)
    rate = if agent.days_infected < model.detection_time
        model.beta_und[agent.pos]
    else
        model.beta_det[agent.pos]
    end

    n_infections = rand(abmrng(model), Poisson(rate))
    n_infections == 0 && return nothing

    neighbors = [a for a in agents_in_position(agent.pos, model) if a.id != agent.id]
    shuffle!(abmrng(model), neighbors)

    for contact in neighbors
        if contact.status == :S
            contact.status = :I
            contact.days_infected = 1
            n_infections -= 1
        elseif contact.status == :R && rand(abmrng(model)) <= model.reinfection_probability
            contact.status = :I
            contact.days_infected = 1
            n_infections -= 1
        end

        n_infections == 0 && return nothing
    end

    return nothing
end

function recover_or_die!(agent, model)
    if agent.status == :I && agent.days_infected >= model.infection_period
        if rand(abmrng(model)) <= model.death_rate
            remove_agent!(agent, model)
        else
            agent.status = :R
            agent.days_infected = 0
        end
    end
    return nothing
end

function safe_step!(model, n::Int = 1)
    for _ in 1:n
        for id in collect(allids(model))
            agent = try
                model[id]
            catch
                nothing
            end
            agent === nothing || sir_agent_step!(agent, model)
        end
        update_quarantine!(model)
    end
    return model
end

infected_count(model) = count(a.status == :I for a in allagents(model))
recovered_count(model) = count(a.status == :R for a in allagents(model))
susceptible_count(model) = count(a.status == :S for a in allagents(model))
total_count(model) = nagents(model)
locked_city_count(model) = count(identity, model.city_locked_down)

function city_status_count(model, city::Integer, status::Symbol)
    return count(a.status == status for a in agents_in_position(city, model))
end

city_population(model, city::Integer) = length(collect(ids_in_position(city, model)))

function infected_fraction(model; denominator = sum(model.Ns))
    denominator == 0 && return 0.0
    return infected_count(model) / denominator
end

function lock_city!(city::Integer, model)
    model.city_locked_down[city] && return nothing
    model.city_locked_down[city] = true
    model.migration_rates[city, :] .= 0.0
    model.migration_rates[city, city] = 1.0
    return nothing
end

function update_quarantine!(model)
    threshold = model.quarantine_threshold
    threshold === nothing && return nothing

    for city in 1:model.C
        model.city_locked_down[city] && continue
        population = city_population(model, city)
        population == 0 && continue
        infected_share = city_status_count(model, city, :I) / population
        infected_share >= threshold && lock_city!(city, model)
    end

    return nothing
end
