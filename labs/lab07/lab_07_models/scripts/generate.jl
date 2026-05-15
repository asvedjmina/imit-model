using DrWatson
@quickactivate "lab_07_models"
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
  cache: false
---

"""
end

sources = [
    (
        "discrete_event_literate.jl",
        "discrete_event",
        "Discrete-event models: M/M/c and Ross model",
        "Base experiments for queueing and repair systems",
    ),
    (
        "discrete_event_param_literate.jl",
        "discrete_event_param",
        "Discrete-event models: parameter sweeps",
        "Server-count, machine-count and repairer-count studies",
    ),
]

for (src_file, base, title, description) in sources
    src_path = joinpath(scripts_dir, src_file)

    Literate.script(src_path, scripts_dir; name = base * "_clean", credit = false)
    Literate.notebook(src_path, notebooks_dir; name = base, execute = false, credit = false)
    Literate.markdown(
        src_path,
        docs_dir;
        name = base,
        flavor = Literate.QuartoFlavor(),
        credit = false,
    )

    qmd_path = joinpath(docs_dir, base * ".qmd")
    content = read(qmd_path, String)
    write(qmd_path, quarto_header(title, description) * content)
end

println("Literate generation completed.")
