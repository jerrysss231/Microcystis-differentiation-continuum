# Global genomic differentiation in *Microcystis*

This repository contains analysis code for a global population-genomic study of *Microcystis*. The analyses examine phylogenomic and population structure, genome-wide diversity and differentiation, homologous recombination and gene flow, genome-content variation, geographic and environmental associations, and projected genomic offset under future climate conditions.

The main analyses are organized by biological question rather than by figure number.

## Repository contents

```text
analysis/
├── population_structure/
│   ├── population_structure_ani.R
│   └── lineage_geography.R
├── population_genetics/
│   └── popgenome_statistics.R
├── differentiation/
│   └── differentiation_continuum.R
├── recombination/
│   ├── shared_recombination_gene_flow.R
│   └── recombination_parameters.R
├── landscape_genomics/
│   └── landscape_genomics.R
├── genomic_offset/
│   ├── gradient_forest_offset.R
│   └── offset_geography.R
└── phylogenomics/
    └── tree_concordance.R

data/
└── README.md

METHODS.md
```

`METHODS.md` gives the analysis settings used in the manuscript, including software versions, statistical tests and environmental data sources.

## Analyses

- **Population structure and ANI**: fastBAPS and snapclust clustering, FastANI summaries, lineage assignments and geographic comparisons.
- **Population-genetic statistics**: PopGenome estimates of nucleotide diversity, pairwise FST and Dxy in 50-kb windows with a 12.5-kb step, and Tajima's D and Fu and Li's F in 10-kb windows with a 2.5-kb step. This workflow was adapted from the population-genomic analysis of Stanojković et al. (2024, *Nature Communications* 15:2122).
- **Differentiation continuum**: pairwise FST, Dxy, ANI distance, gene-content distance and gene-flow resistance; standardized composite scores; PCA; matrix-permutation tests.
- **Homologous recombination**: Gubbins-derived within- and between-lineage shared recombination, gene-flow resistance, rho/theta and r/m comparisons.
- **Landscape genomics**: geographic, environmental, ANI, SNP and gene-content distance matrices; Mantel and partial Mantel tests; dbRDA and variation partitioning.
- **Genomic offset**: Gradient Forest models fitted to LFMM2 candidate environment-associated SNPs, followed by projection to mid-century SSP2-4.5 climate conditions and spatial summaries of genomic offset.
- **Tree concordance**: Robinson-Foulds, information-based tree distances and Mantel comparisons among phylogenomic trees.

## Running the scripts

The scripts use paths relative to a project root. Set the `MICROCYSTIS_PROJECT_ROOT` environment variable to the directory containing `data/` and `results/`:

```bash
export MICROCYSTIS_PROJECT_ROOT=/path/to/Microcystis-project
```

On Windows PowerShell:

```powershell
$env:MICROCYSTIS_PROJECT_ROOT="D:\path\to\Microcystis-project"
```

If the variable is not set, the current working directory is used. Expected input files are listed in [`data/README.md`](data/README.md). Output directories are created under `results/` by the individual scripts.

The R analyses use packages including `PopGenome`, `fastbaps`, `adegenet`, `ape`, `phangorn`, `TreeDist`, `vegan`, `geosphere`, `IRanges`, `gradientForest`, `extendedForest`, `terra`, `FSA`, and the tidyverse packages used in the scripts. Upstream genome-processing software and versions are listed in `METHODS.md`.

## Data availability

Genome assemblies are publicly available from NCBI GenBank and RefSeq; accession numbers are reported with the study metadata. Contemporary bioclimatic and solar-radiation layers were obtained from WorldClim v2.1, and UV-B variables from glUV. Future climate projections use WorldClim v2.1 MPI-ESM1-2-HR data for SSP2-4.5, 2041-2060.

Large genome alignments, external databases and climate rasters are not duplicated in this repository. Derived data files distributed with the article or deposited in a public archive can be placed under the paths described in `data/README.md`.

## Citation

Citation information will be added when the article is published.
