using DrWatson
@quickactivate "lab_02_models"
using Literate

# Директории
scripts_dir = @__DIR__
nb_dir      = joinpath(@__DIR__, "..", "notebooks")
docs_dir    = joinpath(@__DIR__, "..", "docs")

mkpath(nb_dir)
mkpath(docs_dir)

# YAML front matter для Quarto-документов
function quarto_header(title, description)
    """
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
    ("lv_ode_literate.jl",
     "lv_ode",
     "Модель Лотки–Вольтерры (хищник–жертва)",
     "Численное решение и визуализация системы ОДУ «хищник–жертва»"),
    ("sir_ode_literate.jl",
     "sir_ode",
     "Модель SIR (эпидемиологическая)",
     "Численное моделирование распространения инфекции методом SIR"),
]

for (src_file, base, title, description) in sources
    src_path = joinpath(scripts_dir, src_file)
    println("="^60)
    println("Источник: $src_file")

    # 1. Чистый скрипт
    Literate.script(src_path, scripts_dir;
        name   = base * "_clean",
        credit = false)
    println("  ✓ scripts/$(base)_clean.jl")

    # 2. Jupyter notebook (без выполнения — выполним отдельно)
    Literate.notebook(src_path, nb_dir;
        name    = base,
        execute = false,
        credit  = false)
    println("  ✓ notebooks/$base.ipynb")

    # 3. Quarto документ
    qmd_path = joinpath(docs_dir, base * ".qmd")
    Literate.markdown(src_path, docs_dir;
        name   = base,
        flavor = Literate.QuartoFlavor(),
        credit = false)
    # Добавляем YAML front matter в начало
    content = read(qmd_path, String)
    write(qmd_path, quarto_header(title, description) * content)
    println("  ✓ docs/$base.qmd")
end

println("="^60)
println("Генерация завершена.")
println("Файлы:")
println("  scripts/  — чистые скрипты Julia")
println("  notebooks/ — Jupyter-ноутбуки (.ipynb)")
println("  docs/      — Quarto-документация (.qmd)")
