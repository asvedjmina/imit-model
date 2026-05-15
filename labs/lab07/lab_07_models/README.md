# lab_07_models

Седьмая лабораторная работа по курсу "Имитационное моделирование".

Проект включает две модели дискретно-событийного моделирования:

- `M/M/c` для многоканальной системы массового обслуживания;
- модель Росса для отказов, резерва и ремонта машин.

Структура проекта:

- `src/DiscreteEventModels.jl` — общий модуль с симуляциями, аналитикой и построением графиков;
- `scripts/` — literate-источники, clean-скрипты и генератор производных форматов;
- `docs/` — Quarto-документация, сгенерированная из literate-кода;
- `notebooks/` — Jupyter notebooks, сгенерированные из literate-кода;
- `data/` — таблицы с результатами экспериментов;
- `plots/` — графики для обеих моделей;
- `test/` — автоматические тесты.

Основные команды:

```bash
julia --project=. scripts/generate.jl
julia --project=. scripts/discrete_event_clean.jl
julia --project=. scripts/discrete_event_param_clean.jl
julia --project=. test/runtests.jl
```
