# Задание 3 — моделирование орбиты по SGP4

Скрипт `simulate_orbit.py` запускает SGP4 по сохранённому TLE NORAD 57172 на 7 суток с шагом 60 секунд. Основные параметры выводятся в терминал, графики сохраняются в `results/`.

## Запуск

```bash
python -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
python simulate_orbit.py
```

Отчёт собран pdfLaTeX из `report.tex`:

```bash
pdflatex -interaction=nonstopmode report.tex
pdflatex -interaction=nonstopmode report.tex
```

Объект NORAD 57172 — UMKA 1 (RS40S). Его TLE сохранён в `data/tle_57172.txt`. Численные результаты относятся к эпохе этого набора; обычный запуск использует сохранённую копию для воспроизводимости.

Чтобы использовать новый TLE, замените содержимое `data/tle_57172.txt` данными спутника с CelesTrak и запустите скрипт снова. Отчёт содержит результаты для прежнего сохранённого набора; после замены TLE его численные значения нужно обновить.

`questions.md` содержит ответы на частые вопросы защиты.
