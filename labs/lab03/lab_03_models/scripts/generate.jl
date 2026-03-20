using DrWatson
@quickactivate "lab_03_models"
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
    ("daisyworld_literate.jl",
     "daisyworld",
     "Модель Daisyworld: агентное моделирование",
     "Реализация и визуализация агентной модели Daisyworld на Julia с помощью Agents.jl"),
    ("daisyworld_param_literate.jl",
     "daisyworld_param",
     "Модель Daisyworld: исследование параметров",
     "Параметрическое исследование модели Daisyworld"),
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
