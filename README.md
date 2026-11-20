# OEL studies

Can an AI analyst given the same orbital mechanics problem as a published paper recover the authors’ conclusion? This repository collects those experiments. Each paper has a standalone workspace with the inputs, OEL version, replay command and numerical comparison. You can inspect the assumptions, run the simulation yourself and compare your result with the paper.

Orbital Engineering Lab (OEL) is the simulation engine used in these studies. Its source, documentation and installation instructions are available in the [public OEL repository](https://github.com/adamcohen8/orbital-engineering-lab).

## Study index

| Paper and workspace | Question | Verified replay |
| --- | --- | --- |
| [Schaub & Alfriend (2001), *J₂ Invariant Relative Orbits for Spacecraft Formations*](studies/schaub-alfriend-2001/) | How much does initializing a formation in mean orbital elements reduce relative phase drift over 45 orbits? | OEL 1.1.0 measured a **305.5×** reduction in the paper’s 48° example, consistent with Figures 3B and 4B. |

## Replay the first study

Install the pinned OEL version using the study’s instructions. The verified setup used Apple Silicon, CPython 3.11 and the signed OEL 1.1.0 core release. The replay runs locally with OEL and Python’s standard library, without a model API key. [Installation and replay instructions](studies/schaub-alfriend-2001/README.md#replay)

With the OEL launcher on your path, download just the study and run it:

```sh
git clone --depth 1 --filter=blob:none --sparse https://github.com/adamcohen8/OEL-studies.git
cd OEL-studies
git sparse-checkout set studies/schaub-alfriend-2001
cd studies/schaub-alfriend-2001
python3 reproduce.py
```

The script prints the signed phase drifts and their reduction factor, then writes them with a pass flag for each case to `outputs/replay-*/comparison.json`. The 30° argument-of-perigee case reproduces the paper’s example. The 75° case tests whether its conclusion holds when that one input changes. Both use a 48° inclination. [Measurement method and reference values](studies/schaub-alfriend-2001/README.md#measurement-and-comparison)

## How we check a reproduction

1. **Define the claim.** Choose a specific figure, table or conclusion. Record the physical assumptions, units, reference values and acceptable numerical differences. A value estimated from a plotted curve carries less precision than a tabulated value.
2. **Check OEL’s capabilities.** Identify which steps OEL supports, including how the initial conditions are prepared. If an essential capability is missing, pause the study to decide whether it warrants a general OEL feature or a study-specific bridge.
3. **Observe the analyst.** Give GPT-6.2 Sol the problem and withhold the reference answer. Document what the model chose, what the host supplied and what OEL computed. In the first study, for example, the host prepared the starting states and Sol chose simulation settings and measured the drift.
4. **Verify the public workspace.** Replay the simulation from the paper’s folder using the pinned engine. Compare the measured result with the reference and report the scope of agreement and any discrepancies. When a study leads to an OEL improvement, its workspace pins the exact tested development version and immutable revision.
