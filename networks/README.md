# Empirical antagonistic networks

This directory contains the **122 empirical bipartite antagonistic networks**
used in the article:

> Gallarta-Sáenz, P., Lamata-Otín, S., Gómez-Gardeñes, J. & de Andreazzi, C. S.
> *Network structure and species connectivity drive asymmetric fitness
> outcomes in coevolving antagonistic communities.*

The networks are used to evaluate the trait-based coevolutionary model
described in the article (Eqs. 1–6 of the main text) under realistic
ecological conditions, complementing the results obtained on the synthetic
nested network of Fig. 1a.

## Data sources

The networks were compiled from the following sources:

- **Interaction Web Database** (IWDB) — www.nceas.ucsb.edu/interactionweb/index.html
- **Web of Life** — www.web-of-life.es
- Pires & Guimarães (2013)
- Flores et al. (2011)
- Bellay et al. (2013)
- Andreazzi et al. (2020)

The networks span a broad diversity of antagonistic interaction types (e.g.
predator-prey, host-parasite, herbivore-plant), community sizes, and
topological architectures, including both nested and modular networks.

## Data structure

Each network is represented as a **bipartite adjacency matrix** `f_ij`, where
rows correspond to exploiter species and columns to victim species, with
`f_ij = 1` if exploiter species *i* interacts with victim species *j*, and
`f_ij = 0` otherwise.

```
network_XXX/
├── adjacency_matrix.csv    # Binary interaction matrix (exploiters × victims)
├── species_exploiters.csv  # Numeric identifiers of exploiter species (0..M_E-1)
├── species_victims.csv     # Numeric identifiers of victim species (0..M_V-1)
└── metadata.json           # Original network name, interaction type, connectance
```

Networks are numbered `network_001` … `network_122`. Each `metadata.json`
also records the `original_id` (`AN_1` … `AN_122`) used internally during
data processing, and the `original_name` of the source network, so results
can be traced back to the compiled source dataset. Original species names
are not available for these networks, so exploiter/victim species are
identified by numeric index only.

### `network_000`: synthetic illustrative network

`network_000` is a small (5 exploiters x 5 victims), perfectly nested
synthetic network used to illustrate the coevolutionary model (Fig. 1a and
the temporal-dynamics panels of Fig. 2/3), **not** part of the 122-network
empirical dataset above: it is not included in `network_descriptors.csv` or
in any of the PCA/statistical analyses run over `network_001..network_122`.

## Computed structural descriptors

For each network, the eight structural descriptors used in the principal
component analysis (PCA) of Section II.C / III.C of the article were
computed:

| Symbol | Descriptor | Column in `network_descriptors.csv` |
|---|---|---|
| `M_E` | Number of exploiter species | `M_E` |
| `M_V` | Number of victim species | `M_V` |
| `E/V` | Exploiter-to-victim ratio | `E/V` |
| `⟨k_E⟩` | Mean exploiter degree | `mean_k_E` |
| `⟨k_V⟩` | Mean victim degree | `mean_k_V` |
| `ρ_L` | Connectance | `rho_L` |
| `Q` | Modularity | `Q` |
| `NODF` | Nestedness (Almeida-Neto et al. 2008) | `NODF` |

These descriptors are available in `network_descriptors.csv`, with one row
per network matching the identifiers used in `metadata.json` (`id` /
`original_id` columns).

## Simulation code

The code used to simulate the coevolutionary dynamics (Eqs. 1–4) and compute
species fitness (Eqs. 5–6) on these networks is available at:

> *(fill in the link/DOI of the code repository, e.g. Zenodo/GitHub)*

## Citation

If you use these networks, please cite the original article and the primary
source of each individual network, as well as:

> Almeida-Neto, M., Guimarães, P., Guimarães, P. R., Loyola, R. D. & Ulrich, W.
> (2008). A consistent metric for nestedness analysis in ecological systems:
> reconciling concept and measurement. *Oikos*, 117, 1227–1239.

## License

*(specify the data redistribution license, taking into account the terms of
use of IWDB and Web of Life)*
