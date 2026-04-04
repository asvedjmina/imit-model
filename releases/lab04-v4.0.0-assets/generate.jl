using DrWatson
@quickactivate "lab_04_models"
using Literate

scripts_dir = @__DIR__
notebooks_dir = joinpath(@__DIR__, "..", "notebooks")
docs_dir = joinpath(@__DIR__, "..", "docs")

mkpath(notebooks_dir)
mkpath(docs_dir)

function quarto_header(title, description)
    return """
---
title: "$(title)"
description: "$(description)"
lang: ru-RU
jupyter: julia-1.12
format:
  html:
    toc: true
    number-sections: true
    code-fold: false
execute:
  cache: true
---

"""
end

sources = [
    (
        "sir_literate.jl",
        "sir",
        "SIR model on a graph of cities",
        "Base experiments for the graph-based SIR model",
    ),
    (
        "sir_param_literate.jl",
        "sir_param",
        "SIR model: parameter scenarios",
        "Scenario analysis for the graph-based SIR model",
    ),
    (
        "sir_threshold_literate.jl",
        "sir_threshold",
        "SIR model: threshold study",
        "Threshold study for the graph-based SIR model",
    ),
    (
        "sir_heterogeneity_literate.jl",
        "sir_heterogeneity",
        "SIR model: heterogeneity effect",
        "Heterogeneity study for the graph-based SIR model",
    ),
    (
        "sir_quarantine_literate.jl",
        "sir_quarantine",
        "SIR model: quarantine policy",
        "Quarantine policy analysis for the graph-based SIR model",
    ),
    (
        "sir_optimize_literate.jl",
        "sir_optimize",
        "SIR model: constrained optimization",
        "Constrained optimization for the graph-based SIR model",
    ),
]

for (src_file, base, title, description) in sources
    src_path = joinpath(scripts_dir, src_file)

    Literate.script(src_path, scripts_dir; name = base * "_clean", credit = false)
    Literate.notebook(src_path, notebooks_dir; name = base, execute = true, credit = false)
    Literate.markdown(src_path, docs_dir; name = base, flavor = Literate.QuartoFlavor(), credit = false)

    qmd_path = joinpath(docs_dir, base * ".qmd")
    content = read(qmd_path, String)
    write(qmd_path, quarto_header(title, description) * content)
end

println("Literate generation completed.")
