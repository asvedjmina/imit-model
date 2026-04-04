# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [4.0.0] - 2026-04-04

### Added
- Проект DrWatson `lab_04_models` для лабораторной работы №4
- Агентная SIR-модель на Agents.jl для трёх городов
- Базовый сценарий, сканирование `beta`, исследование миграции и сценарный анализ
- Дополнительные эксперименты: порог эпидемии, гетерогенность, карантин и оптимизация
- Шесть literate-скриптов с генерацией `.jl`, `.ipynb`, `.qmd`
- Выполненные Jupyter notebooks для всех сценариев лабораторной работы №4
- Отчёт по лабораторной работе №4 (Quarto, PDF, DOCX)
- Презентация по лабораторной работе №4 (Beamer, reveal.js)

### Changed
- Обновлены библиография и LaTeX-преамбула отчёта для корректной сборки и ссылок
- Оптимизационный сценарий приведён к многокритериальной постановке через `BlackBoxOptim`

## [3.0.0] - 2026-03-20

### Added
- Проект DrWatson `lab_03_models` для лабораторной работы №3
- Агентная модель Daisyworld на Agents.jl
- Литературные Julia-скрипты `daisyworld_literate.jl` и `daisyworld_param_literate.jl`
- Генерация производных форматов через Literate.jl (`.jl`, `.ipynb`, `.qmd`)
- Выполненные Jupyter notebooks для базового и параметрического сценариев
- Визуализации Daisyworld: heatmap, динамика численности, luminosity ramp, анимация
- Отчёт по лабораторной работе №3 (Quarto, PDF, DOCX)
- Презентация по лабораторной работе №3 (Beamer, reveal.js)

### Changed
- Обновлён bib-файл: добавлены источники по агентному моделированию и Daisyworld

## [2.0.0] - 2026-03-07

### Added
- Проект DrWatson `lab_02_models` для лабораторной работы №2
- Литературный Julia-скрипт `sir_ode_literate.jl` (модель SIR)
- Литературный Julia-скрипт `lv_ode_literate.jl` (модель Лотки-Вольтерры)
- Генерация производных форматов через Literate.jl (.jl, .ipynb, .qmd)
- Выполненные Jupyter notebooks для обеих моделей
- Анализ чувствительности к параметрам (SIR и LV)
- Отчёт по лабораторной работе №2 (Quarto, PDF, DOCX)
- Презентация по лабораторной работе №2 (Beamer, reveal.js)

### Changed
- Обновлён bib-файл: добавлены ссылки на SIR, LV, Literate.jl, DrWatson

## [1.0.0] - 2026-02-12

### Added
- Структура курса имитационного моделирования
- Лабораторная работа №1
