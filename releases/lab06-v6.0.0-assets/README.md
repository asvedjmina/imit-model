# lab_06_models

Шестая лабораторная работа по курсу "Имитационное моделирование".

Структура проекта:

- `src/SIRPetri.jl` — модуль с реализацией SIR-модели в терминах сети Петри;
- `scripts/` — literate-источники, clean-скрипты и генератор форматов;
- `docs/` — Quarto-документация, сгенерированная из literate-кода;
- `notebooks/` — Jupyter notebook, сгенерированные из literate-кода;
- `data/` — таблицы с результатами экспериментов;
- `plots/` — графики и GIF-анимация;
- `test/` — минимальные автоматические тесты.

Основные команды:

```bash
julia --project=. scripts/generate.jl
julia --project=. scripts/sirpetri_clean.jl
julia --project=. scripts/sirpetri_param_clean.jl
julia --project=. test/runtests.jl
```
