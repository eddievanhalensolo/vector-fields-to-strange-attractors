# Dynamische Systeme & Chaos

Eine interaktive Lernreise von Differentialgleichungen bis zur Chaos-Theorie.

## Struktur

```
dynamical-systems-chaos/
├── notebooks/
│   ├── 01_dgl_als_vektorfeld.ipynb      # Grundlagen: DGL als Vektorfeld
│   ├── 02_eindimensionale_systeme.ipynb  # 1D: Fixpunkte, Stabilität, Potentiale
│   ├── 03_zweidimensionale_systeme.ipynb # 2D: Phasenportraits, Klassifikation
│   ├── 04_nichtlineare_systeme.ipynb     # Linearisierung, Grenzzyklen
│   ├── 05_bifurkationen.ipynb            # Qualitative Übergänge
│   ├── 06_chaos.ipynb                    # Lorenz, Lyapunov, Schmetterlingseffekt
│   └── 07_seltsame_attraktoren.ipynb     # Fraktale Dimension
├── utils/
│   ├── plotting.py                        # Visualisierungsfunktionen
│   └── integration.py                     # Numerische Werkzeuge
└── requirements.txt
```

## Die Lernreise

```
[1] DGL als Vektorfeld
         ↓
[2] 1D-Systeme (kein Chaos möglich!)
         ↓
[3] 2D-Systeme (immer noch kein Chaos)
         ↓
[4] Nichtlinearität & Grenzzyklen
         ↓
[5] Bifurkationen
         ↓
[6] Chaos (≥3D oder diskret)
         ↓
[7] Seltsame Attraktoren
```

## Setup

```bash
cd dynamical-systems-chaos
pip install -r requirements.txt
jupyter notebook
```

## Kernkonzepte

- **DGL = Vektorfeld**: Jede DGL definiert einen Fluss im Phasenraum
- **Phasenraum-Denken**: Qualitative Analyse ohne explizite Lösungen
- **Linearisierung**: Lokal verhält sich alles wie linear (Hartman-Grobman)
- **Bifurkationen**: Qualitative Übergänge bei Parameteränderung
- **Chaos**: Determinismus ohne Vorhersagbarkeit (λ > 0)
- **Fraktale Dimension**: Die Geometrie seltsamer Attraktoren

## Literatur

- **Strogatz**: "Nonlinear Dynamics and Chaos" (das Standardwerk)
- **Gleick**: "Chaos: Making a New Science" (populärwissenschaftlich)
- **Lorenz**: "The Essence of Chaos" (vom Entdecker selbst)
