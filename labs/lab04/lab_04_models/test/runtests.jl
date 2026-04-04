using DrWatson
using Test

@quickactivate "lab_04_models"

include(srcdir("sir_model.jl"))

println("Starting tests")

@testset "lab_04_models tests" begin
    migration_matrix = create_migration_matrix(3, 0.2)
    @test all(isapprox(sum(migration_matrix[i, :]), 1.0; atol = 1e-8) for i in 1:3)
    @test all(isapprox(migration_matrix[i, i], 0.8; atol = 1e-8) for i in 1:3)

    model = initialize_sir(
        Ns = [10, 10, 10],
        beta_und = [0.0, 0.0, 0.0],
        beta_det = [0.0, 0.0, 0.0],
        infection_period = 3,
        detection_time = 1,
        death_rate = 0.0,
        reinfection_probability = 0.0,
        Is = [0, 0, 2],
        migration_rates = create_migration_matrix(3, 0.0),
        seed = 7,
        n_steps = 5,
    )

    @test total_count(model) == 30
    @test infected_count(model) == 2
    @test susceptible_count(model) == 28
    @test recovered_count(model) == 0

    safe_step!(model, 3)
    @test infected_count(model) == 0
    @test recovered_count(model) == 2
    @test total_count(model) == 30

    quarantine_model = initialize_sir(
        Ns = [10, 10, 10],
        beta_und = [0.0, 0.0, 0.0],
        beta_det = [0.0, 0.0, 0.0],
        infection_period = 5,
        detection_time = 1,
        death_rate = 0.0,
        reinfection_probability = 0.0,
        Is = [1, 0, 0],
        migration_rates = create_migration_matrix(3, 0.2),
        quarantine_threshold = 0.01,
        seed = 11,
        n_steps = 3,
    )

    safe_step!(quarantine_model, 1)
    @test quarantine_model.city_locked_down[1]
    @test quarantine_model.migration_rates[1, 1] == 1.0
    @test sum(quarantine_model.migration_rates[1, :]) == 1.0
end
