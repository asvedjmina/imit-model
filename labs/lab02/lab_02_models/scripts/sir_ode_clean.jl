using DrWatson
@quickactivate "lab_02_models"
using DifferentialEquations
using SimpleDiffEq
using Tables
using DataFrames
using StatsPlots
using LaTeXStrings
using Plots
using BenchmarkTools

script_name = splitext(basename(PROGRAM_FILE))[1]
mkpath(plotsdir(script_name))
mkpath(datadir(script_name))

function sir_ode!(du, u, p, t)
    (S, I, R) = u
    (β, c, γ) = p
    N = S + I + R
    @inbounds begin
        du[1] = -β * c * I / N * S          # dS/dt
        du[2] =  β * c * I / N * S - γ * I  # dI/dt
        du[3] =  γ * I                       # dR/dt
    end
    nothing
end

δt    = 0.1
tmax  = 40.0
tspan = (0.0, tmax)
u0    = [990.0, 10.0, 0.0]   # [S₀, I₀, R₀]
p     = [0.05, 10.0, 0.25]   # [β, c, γ]

R0 = (p[2] * p[1]) / p[3]   # R₀ = (c·β) / γ

println("="^60)
println("Параметры модели SIR")
println("="^60)
println("β = $(p[1]),  c = $(p[2]),  γ = $(p[3])")
println("R₀ = c·β/γ = $(round(R0, digits=3))")
println("Средняя продолжительность болезни = $(round(1/p[3], digits=2)) дней")
println("Начальные условия: S₀=$(u0[1]), I₀=$(u0[2]), R₀=$(u0[3])")

prob_ode = ODEProblem(sir_ode!, u0, tspan, p)
sol_ode  = solve(prob_ode, dt = δt)

df_ode = DataFrame(Tables.table(sol_ode'))
rename!(df_ode, ["S", "I", "R"])
df_ode[!, :t] = sol_ode.t
df_ode[!, :N] = df_ode.S + df_ode.I + df_ode.R   # контроль сохранения N

plt1 = @df df_ode plot(:t, [:S :I :R],
    label     = [L"S(t)" L"I(t)" L"R(t)"],
    xlabel    = "Время, дни",
    ylabel    = "Количество людей",
    title     = "Модель SIR: Динамика эпидемии",
    linewidth = 2,
    legend    = :right,
    grid      = true,
    size      = (800, 500))
annotate!(plt1, maximum(df_ode.t) * 0.7, maximum(df_ode.N) * 0.8,
    text("β=$(p[1]),  c=$(p[2]),  γ=$(p[3])\nR₀=$(round(R0, digits=2))", 8, :left))

savefig(plt1, plotsdir(script_name, "sir_main.png"))

peak_idx   = argmax(df_ode.I)
peak_time  = df_ode.t[peak_idx]
peak_value = df_ode.I[peak_idx]

plt2 = @df df_ode plot(:t, :I,
    label     = L"I(t)",
    xlabel    = "Время, дни",
    ylabel    = "Количество инфицированных",
    title     = "Динамика числа заражённых",
    color     = :red,
    linewidth = 2,
    fill      = (0, 0.3, :red),
    grid      = true,
    size      = (800, 400))
vline!(plt2, [peak_time], color=:black, linestyle=:dash, label=false, linewidth=1)
annotate!(plt2, peak_time, peak_value * 1.05,
    text("Пик: $(round(peak_value, digits=1)) чел.\nна $(round(peak_time, digits=1)) день",
    8, :top))

savefig(plt2, plotsdir(script_name, "sir_infected.png"))

plt3 = @df df_ode plot(:t, :I,
    label     = L"I(t)",
    xlabel    = "Время, дни",
    ylabel    = "Инфицированные (лог. масштаб)",
    title     = "Экспоненциальный рост (лог. шкала)",
    yscale    = :log10,
    color     = :red,
    linewidth = 2,
    grid      = true,
    size      = (800, 400))

savefig(plt3, plotsdir(script_name, "sir_log_scale.png"))

plt4 = @df df_ode plot(:t, [:S :I :R] ./ df_ode.N .* 100,
    label     = [L"S(t)/N" L"I(t)/N" L"R(t)/N"],
    xlabel    = "Время, дни",
    ylabel    = "Доля популяции, %",
    title     = "Динамика эпидемии (в процентах)",
    linewidth = 2,
    legend    = :right,
    grid      = true,
    size      = (800, 500))
if R0 > 1
    herd_immunity_threshold = (1 - 1/R0) * 100
    hline!(plt4, [herd_immunity_threshold], color=:purple, linestyle=:dash,
        label="Порог коллективного иммунитета ($(round(herd_immunity_threshold, digits=1))%)",
        linewidth=1.5)
end

savefig(plt4, plotsdir(script_name, "sir_percentages.png"))

plt5 = plot(df_ode.S, df_ode.I,
    label     = "Фазовая траектория",
    xlabel    = L"S(t)",
    ylabel    = L"I(t)",
    title     = "Фазовый портрет SIR модели",
    color     = :blue,
    linewidth = 2,
    grid      = true,
    size      = (800, 500),
    legend    = :topright)
for i in 1:50:length(df_ode.S)-1
    plot!(plt5, [df_ode.S[i], df_ode.S[i+1]], [df_ode.I[i], df_ode.I[i+1]],
        arrow=:closed, color=:blue, alpha=0.5, label=false)
end

savefig(plt5, plotsdir(script_name, "sir_phase_portrait.png"))

df_ode[!, :Re] = R0 .* df_ode.S ./ df_ode.N

plt6 = @df df_ode plot(:t, :Re,
    label     = L"R_e(t)",
    xlabel    = "Время, дни",
    ylabel    = L"R_e",
    title     = "Динамика эффективного репродуктивного числа",
    color     = :green,
    linewidth = 2,
    grid      = true,
    size      = (800, 400))
hline!(plt6, [1.0], color=:red, linestyle=:dash,
    label="Порог эпидемии (Rₑ = 1)", linewidth=1.5)

cross_idx = findfirst(x -> x < 1, df_ode.Re)
if !isnothing(cross_idx) && cross_idx > 1
    cross_time = df_ode.t[cross_idx]
    vline!(plt6, [cross_time], color=:black, linestyle=:dash, label=false, linewidth=1)
    annotate!(plt6, cross_time, 1.2,
        text("Rₑ < 1 с $(round(cross_time, digits=1)) дня", 8, :left))
end

savefig(plt6, plotsdir(script_name, "sir_effective_R.png"))

plt7 = plot(layout=(2, 3), size=(1200, 800))
plot!(plt7[1], df_ode.t, df_ode.S, label=L"S(t)", color=1, linewidth=2, title="Восприимчивые")
plot!(plt7[2], df_ode.t, df_ode.I, label=L"I(t)", color=2, linewidth=2, title="Заражённые")
plot!(plt7[3], df_ode.t, df_ode.R, label=L"R(t)", color=3, linewidth=2, title="Выздоровевшие")
plot!(plt7[4], df_ode.t, df_ode.I, label=L"I(t)", color=2, linewidth=2,
    yscale=:log10, title="Лог. масштаб")
plot!(plt7[5], df_ode.S, df_ode.I, label=false, color=4, linewidth=2,
    title="Фазовый портрет", xlabel=L"S", ylabel=L"I")
plot!(plt7[6], df_ode.t, df_ode.Re, label=L"R_e", color=:green, linewidth=2,
    title=L"R_e(t)")
hline!(plt7[6], [1.0], color=:red, linestyle=:dash, label=false)

savefig(plt7, plotsdir(script_name, "sir_panel.png"))

println("\n=== АНАЛИЗ РЕЗУЛЬТАТОВ ===")
println("Общая численность популяции (контроль сохранения): N = $(round(df_ode.N[1], digits=1))")
println("Пиковое число заражённых:    I_max = $(round(peak_value,           digits=1))")
println("Время достижения пика:       t_peak = $(round(peak_time,           digits=1)) дней")
println("Итоговое число переболевших: R(∞)  = $(round(df_ode.R[end],       digits=1))")
println("Доля переболевших:           $(round(df_ode.R[end]/df_ode.N[1]*100, digits=1))%")
if R0 > 1
    println("\nТеоретический анализ:")
    println("  Порог коллективного иммунитета: $(round((1-1/R0)*100, digits=1))%")
    println("  Пик ожидается при S/N = 1/R₀ = $(round(1/R0, digits=3))")
end

param_sets_sir = [
    (label="R₀=0.5",  β=0.0125, c=10.0, γ=0.25),
    (label="R₀=1.0",  β=0.0250, c=10.0, γ=0.25),
    (label="R₀=2.0",  β=0.0500, c=10.0, γ=0.25),
    (label="R₀=3.0",  β=0.0750, c=10.0, γ=0.25),
    (label="R₀=4.0",  β=0.1000, c=10.0, γ=0.25),
]

sweep_results_sir = DataFrame(
    сценарий         = String[],
    R0               = Float64[],
    I_max            = Float64[],
    t_peak           = Float64[],
    R_final          = Float64[],
    охват_процент    = Float64[],
    порог_иммунитета = Float64[],
)

solutions_sir = []
N_base = u0[1] + u0[2] + u0[3]

for ps in param_sets_sir
    p_test  = [ps.β, ps.c, ps.γ]
    R0_test = ps.c * ps.β / ps.γ

    prob_test = ODEProblem(sir_ode!, u0, tspan, p_test)
    sol_test  = solve(prob_test, dt=δt)

    df_test = DataFrame(Tables.table(sol_test'))
    rename!(df_test, ["S", "I", "R"])
    df_test[!, :t]  = sol_test.t
    df_test[!, :Re] = R0_test .* df_test.S ./ N_base

    I_max_test  = maximum(df_test.I)
    t_peak_test = df_test.t[argmax(df_test.I)]
    R_fin_test  = df_test.R[end]
    herd_test   = R0_test > 1 ? (1 - 1/R0_test) * 100 : 0.0

    push!(solutions_sir, (
        t  = df_test.t,
        S  = df_test.S,
        I  = df_test.I,
        R  = df_test.R,
        Re = df_test.Re,
        label = ps.label,
        R0    = R0_test,
    ))
    push!(sweep_results_sir, (
        ps.label,
        round(R0_test,               digits=2),
        round(I_max_test,            digits=1),
        round(t_peak_test,           digits=1),
        round(R_fin_test,            digits=1),
        round(R_fin_test/N_base*100, digits=1),
        round(herd_test,             digits=1),
    ))
end

println("\n" * "="^60)
println("Результаты анализа чувствительности SIR по R₀")
println("="^60)
println(sweep_results_sir)

plt_r1 = plot(xlabel="Время, дни", ylabel="Инфицированные I(t)",
    title="Динамика заражённых при разных R₀",
    grid=true, size=(900, 500), legend=:topright)
for (i, s) in enumerate(solutions_sir)
    plot!(plt_r1, s.t, s.I, label=s.label, linewidth=2, color=i)
end

savefig(plt_r1, plotsdir(script_name, "sir_sweep_infected.png"))

plt_r2 = plot(xlabel="Время, дни", ylabel=L"R_e(t)",
    title="Эффективное репродуктивное число при разных R₀",
    grid=true, size=(900, 500), legend=:topright)
hline!(plt_r2, [1.0], color=:black, linestyle=:dash, linewidth=1.5, label="Rₑ = 1")
for (i, s) in enumerate(solutions_sir)
    plot!(plt_r2, s.t, s.Re, label=s.label, linewidth=2, color=i)
end

savefig(plt_r2, plotsdir(script_name, "sir_sweep_Re.png"))

R0_vals    = sweep_results_sir.R0
I_max_vals = sweep_results_sir.I_max
cover_vals = sweep_results_sir.охват_процент

plt_r3 = plot(layout=(1, 2), size=(1000, 420))
plot!(plt_r3[1], R0_vals, I_max_vals,
    xlabel=L"R_0", ylabel="Пиковое число заражённых",
    title=L"I_{max} \; vs \; R_0",
    marker=:circle, markersize=7, linewidth=2, color=:red,
    legend=false, grid=true)
vline!(plt_r3[1], [1.0], color=:black, linestyle=:dash, label=false)

plot!(plt_r3[2], R0_vals, cover_vals,
    xlabel=L"R_0", ylabel="Доля переболевших, %",
    title=L"R(\infty)/N \; vs \; R_0",
    marker=:circle, markersize=7, linewidth=2, color=:purple,
    legend=false, grid=true)
vline!(plt_r3[2], [1.0], color=:black, linestyle=:dash, label=false)

savefig(plt_r3, plotsdir(script_name, "sir_sweep_outcomes.png"))

plt_r4 = plot(layout=(2, 2), size=(1200, 800))
for (i, s) in enumerate(solutions_sir)
    plot!(plt_r4[1], s.t, s.I,
        label=s.label, color=i, linewidth=1.8,
        title="I(t)", xlabel="дни", ylabel="I", legend=:topright, grid=true)
    plot!(plt_r4[2], s.t, s.Re,
        label=false, color=i, linewidth=1.8,
        title=L"R_e(t)", xlabel="дни", ylabel=L"R_e", grid=true)
    plot!(plt_r4[3], s.S, s.I,
        label=false, color=i, linewidth=1.8,
        title="Фазовый портрет", xlabel="S", ylabel="I", grid=true)
end
hline!(plt_r4[2], [1.0], color=:black, linestyle=:dash, label=false)
plot!(plt_r4[4], R0_vals, cover_vals,
    xlabel=L"R_0", ylabel="Охват, %", title=L"R(\infty)/N \; vs \; R_0",
    marker=:circle, markersize=7, linewidth=2, color=:purple,
    legend=false, grid=true)
vline!(plt_r4[4], [1.0], color=:black, linestyle=:dash, label=false)

savefig(plt_r4, plotsdir(script_name, "sir_sweep_panel.png"))

@benchmark solve(prob_ode, dt = δt)
