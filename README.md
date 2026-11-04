# From Vector Fields to Strange Attractors

[![CI](https://github.com/eddievanhalensolo/vector-fields-to-strange-attractors/actions/workflows/ci.yml/badge.svg)](https://github.com/eddievanhalensolo/vector-fields-to-strange-attractors/actions/workflows/ci.yml)
[![Lizenz: MIT](https://img.shields.io/badge/Lizenz-MIT-blue.svg)](LICENSE)

Ein Kurs über dynamische Systeme und Chaos in sieben Jupyter-Notebooks. Er beginnt mit Differentialgleichungen als Vektorfeldern und endet bei seltsamen Attraktoren. Jedes Modul enthält Erklärungen, ausführbaren Code, Abbildungen und Übungen.

![Der Lorenz-Attraktor](assets/lorenz.png)

## Module

Die Notebooks enthalten ihre Ergebnisse und lassen sich direkt auf GitHub lesen. Am besten arbeitest du sie der Reihe nach durch.

| # | Modul | Inhalt |
|---|-------|--------|
| 1 | [DGL als Vektorfeld](notebooks/01_dgl_als_vektorfeld.ipynb) | Phasenraum, Eindeutigkeit, Fluss |
| 2 | [Eindimensionale Systeme](notebooks/02_eindimensionale_systeme.ipynb) | Fixpunkte, Stabilität, Potentiale |
| 3 | [Zweidimensionale Systeme](notebooks/03_zweidimensionale_systeme.ipynb) | Phasenportraits, Klassifikation, Nullklinen |
| 4 | [Nichtlineare Systeme](notebooks/04_nichtlineare_systeme.ipynb) | Linearisierung, Grenzzyklen, Satz von Poincaré-Bendixson |
| 5 | [Bifurkationen](notebooks/05_bifurkationen.ipynb) | Sattelknoten, Heugabel, Hopf, Periodenverdopplung |
| 6 | [Chaos](notebooks/06_chaos.ipynb) | Lorenz-System, Lyapunov-Exponenten, logistische Abbildung |
| 7 | [Seltsame Attraktoren](notebooks/07_seltsame_attraktoren.ipynb) | Fraktale Dimension, Hénon, Rössler |

In einer oder zwei Dimensionen gibt es bei kontinuierlichen Systemen kein Chaos. Dafür braucht es mindestens drei Dimensionen oder eine diskrete Abbildung. Der Kurs baut deshalb von 1D über 2D bis zu den chaotischen Systemen in den letzten beiden Modulen auf.

## Installation

```bash
git clone https://github.com/eddievanhalensolo/vector-fields-to-strange-attractors.git
cd vector-fields-to-strange-attractors
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
jupyter notebook notebooks/
```

Getestet mit Python 3.10 und 3.14.

## Aufbau des Repositorys

```
vector-fields-to-strange-attractors/
├── notebooks/         # Die sieben Module
├── utils/             # Optionale Hilfsfunktionen
│   ├── plotting.py        # Abbildungen
│   └── integration.py     # Numerische Werkzeuge
├── tests/             # Tests für utils/
└── requirements.txt
```

Die Notebooks enthalten ihren gesamten Code selbst und brauchen `utils/` nicht. Die Hilfsfunktionen kannst du zum Beispiel für die Übungen verwenden. Vom Hauptverzeichnis aus:

```python
from utils import plot_phase_portrait, lyapunov_exponents_flow
```

Die Tests vergleichen die numerischen Funktionen mit bekannten Ergebnissen:

```bash
pip install pytest
python -m pytest tests
```

## Literatur

- **Strogatz**: "Nonlinear Dynamics and Chaos" (das Standardlehrbuch)
- **Gleick**: "Chaos: Making a New Science" (populärwissenschaftlich)
- **Lorenz**: "The Essence of Chaos" (vom Entdecker des Lorenz-Systems)
- **Grassberger & Procaccia**: "Characterization of Strange Attractors" (Phys. Rev. Lett. 50, 1983; zur Korrelationsdimension in Modul 7)

## Lizenz

[MIT](LICENSE)
