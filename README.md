# Coevolution in antagonistic networks

Simulation code and analysis notebook for:

> Gallarta-Sáenz, P., Lamata-Otín, S., Gómez-Gardeñes, J. & de Andreazzi, C. S.
> *Network structure and species connectivity drive asymmetric fitness
> outcomes in coevolving antagonistic communities.*

This repository reproduces every figure in the article: the trait-based
coevolutionary dynamics (Eqs. 1–4) and species fitness (Eqs. 5–6), simulated
on a synthetic illustrative network (Figs. 2–3) and on 122 empirical
antagonistic networks (Figs. 4–7 and the associated statistical analyses).

## Repository layout

```
src/
  coevolution_antagonistic_networks.c   # simulation code (C)
  run_coevolution_networks.sh           # compiles the C code and runs it over a range of networks
networks/
  network_000/                          # synthetic nested 5x5 network (Figs. 2-3)
  network_001 ... network_122/          # empirical networks (Figs. 4-7), see networks/README.md
  network_descriptors.csv               # per-network structural descriptors used in the PCA
figures.ipynb                           # produces every figure and statistical table from the simulation output
figures/                                # figures and result tables as published (svg/pdf/png/csv)
dataframes/                             # tidy per-species dataframe built by the notebook (gitignored, regenerated on run)
files/                                  # raw simulation output consumed by the notebook (gitignored, see below)
```

`figures/` and `figures.ipynb`'s cell outputs are committed as-is so the
results are visible without running anything. `files/` (raw simulation
output) and `dataframes/df_networks.csv` (built from it) are **not**
committed — see [Reproducing the results](#reproducing-the-results).

## Requirements

- A C compiler (`gcc`; the run script also works with `clang` aliased to
  `gcc`), `bash`, `make`-free.
- Python 3.10+ with the packages in `requirements.txt`:

  ```bash
  pip install -r requirements.txt
  ```

## Reproducing the results

The simulation output is not committed to this repository: at the article's
defaults (100 simulations x 17x17 = 289 `(xi_E, xi_V)` combinations x 122
networks, with the full time series kept for the illustrative network) it is
well over 200 GB, far past what git/GitHub can host. To regenerate it:

### 1. Simulate

```bash
# Empirical networks (Figs. 4, 5, 7 and their statistical analyses):
# final state only, over the full (xi_E, xi_V) grid, for all 122 networks
src/run_coevolution_networks.sh -i 1 -e 122 -t 0

# Synthetic network_000 (Figs. 2-3):
src/run_coevolution_networks.sh -i 0 -e 0 -t 0   # heatmaps (Fig. 3a-d)
src/run_coevolution_networks.sh -i 0 -e 0 -t 1   # z(t) trajectories (Figs. 2, 3e-f)
```

`run_coevolution_networks.sh` compiles `src/coevolution_antagonistic_networks.c`
and runs the resulting binary in parallel over the requested range of
network ids, writing to `files/networks/` (and `files/networks/temporal/`
for `-t 1`). Run `src/run_coevolution_networks.sh -h` for all flags
(`-i`/`-e` network id range, `-t` final-state vs. full time series, `-j` max
parallel jobs).

The full `(xi_E, xi_V)` sweep with `-t 1` (flag_temporal=1) is 289 x N_Sim x
TMAX/h rows per species — hundreds of millions of rows at the defaults. The
Figs. 2/3 temporal panels only ever need 4-6 specific corners
(`xi_E, xi_V` in `{0.1, 0.3, 0.7, 0.9}`), so it is much cheaper to narrow the
sweep before compiling for that run only: in the `CHANGE WHEN NECESSARY`
block of `src/coevolution_antagonistic_networks.c`, set
`xi_d_V_delta = xi_d_E_delta = 0.2` (grid `{0.1, 0.3, 0.5, 0.7, 0.9}`), then
re-run `-t 1` for `network_000`. `run_coevolution_networks.sh` recompiles the
binary on every invocation, so the edit takes effect on the next run.

Simulation parameters (`N_Sim`, `alpha`, `TMAX`, integration step `h`, and
the `(xi_E, xi_V)` grid) are set in that same `CHANGE WHEN NECESSARY` block
at the top of `main()` and require recompiling to change — there is no
runtime flag for them.

### 2. Run the notebook

```bash
jupyter nbconvert --to notebook --execute --inplace figures.ipynb
```

or open `figures.ipynb` in Jupyter and run all cells. It reads the CSVs
under `files/networks/`, builds `dataframes/df_networks.csv`, and writes
every figure and result table to `figures/`.

## Notebook structure

`figures.ipynb` is organized in two parts:

1. **Figures 2 & 3** — `z(t)` trajectories and fitness heatmaps on the
   synthetic `network_000`.
2. **Figures 4, 5, 7 and statistical analyses** — built from the 122
   empirical networks: fitness-vs-degree (Fig. 4, with a cluster-robust
   regression block quantifying the relationship per `(xi_E, xi_V)` corner
   and gremio), PCA of network structural descriptors (Figs. 5, 7), and a
   PCA-fitness statistical block (correlations, R² map over the full grid,
   permutation test, robustness check).

## Data

`networks/` contains the 122 empirical bipartite networks and their
structural descriptors; see `networks/README.md` for sources, format, and
citation requirements for the underlying data.
