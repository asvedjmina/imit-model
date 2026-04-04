# lab_04_models

Проект для лабораторной работы N4 по курсу "Имитационное моделирование".

Содержимое:

- `src/sir_model.jl` - агентная SIR-модель на графе городов.
- `src/sir_analysis.jl` - функции для базового эксперимента и анализа параметров.
- `scripts/` - исполняемые сценарии и литературные исходники.
- `docs/` - Quarto-документация, сгенерированная из literate-кода.
- `notebooks/` - Jupyter notebooks, сгенерированные из literate-кода.
- `plots/` - графики экспериментов.
- `data/` - CSV и JLD2 с результатами.
- `test/` - минимальные тесты модели.

Основные команды:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. scripts/generate.jl
julia --project=. scripts/sir_run_basic.jl
julia --project=. scripts/sir_scan_beta.jl
julia --project=. scripts/sir_migration_effect.jl
julia --project=. scripts/sir_visualize_dynamics.jl
julia --project=. test/runtests.jl
```
