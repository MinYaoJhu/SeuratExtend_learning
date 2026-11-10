# SeuratExtend: study notes <img src="assets/logo.svg" align="right" width="130" alt="SeuratExtend study notes badge"/>

**My worked tutorials for single-cell RNA-seq analysis with Seurat and SeuratExtend**

[![SeuratExtend](https://img.shields.io/badge/SeuratExtend-huayc09.github.io-0F3B3A)](https://huayc09.github.io/SeuratExtend/)
![Seurat](https://img.shields.io/badge/Seurat-v5-blue)
![R](https://img.shields.io/badge/R-4.5-276DC3?logo=r&logoColor=white)

In autumn 2025 I worked through the standard Seurat workflow and then the [SeuratExtend](https://github.com/huayc09/SeuratExtend) tutorials by Yichao Hua, which add publication-quality plotting, functional analysis and an R interface to Python trajectory tools (scVelo, Palantir, CellRank). Each notebook (`.Rmd`) follows a published tutorial with my notes and edits; the rendered `.md` files can be read directly on GitHub.

## 📚 Notebooks

| # | Notebook | Based on | Topics |
|---|---|---|---|
| 0 | [Seurat PBMC 3k tutorial](0_seurat5_pbmc3k_tutorial.md) | [Seurat guided clustering tutorial](https://satijalab.org/seurat/articles/pbmc3k_tutorial) (Satija Lab) | QC, normalization, variable features, PCA, clustering, UMAP, marker genes; about half of this notebook is my own notes and code |
| 1 | [Quick start-up guide](1_Quick_Start-Up_Guide.md) | SeuratExtend "Quick start" vignette | Cluster and composition plots, marker heatmaps and dot plots, volcano plots, GSEA, SCENIC, and a first look at Palantir and scVelo |
| 2 | [Enhanced visualization](2_Enhanced_Visualization.md) | SeuratExtend "Visualization" vignette | `DimPlot2`, feature plots, heatmaps, violin and box plots, colour palettes (`color_pro`) |
| 3 | [Trajectory and pseudotime analysis](3_Trajectory_and_Pseudotime_Analysis.md) | SeuratExtend "Trajectory" vignette | RNA velocity (scVelo), pseudotime and cell fate (Palantir), MAGIC imputation, Slingshot trajectories and gene-trend plots, Seurat ↔ AnnData conversion |

**Notebook 3** was rendered in November 2026 (R 4.5.3, Seurat 5.5.1, SeuratExtend 1.3.2, with SeuratExtend's Python environment: scVelo 0.3.0, Palantir 1.3.0, CellRank 2.0.0). As in the original SeuratExtend vignette, a few chunks are shown as code only: interactive cell selection, the CellRank plots and a final AnnData plot. Calling CellRank from R Markdown hangs on this Linux system (a known issue when Python's multithreaded libraries run inside R), so those two plots are not included. Notebooks 0–2 were rendered in 2025 with R 4.5.1, Seurat 5.3.0 and SeuratExtend 1.2.5.

## 🧠 Skills practiced

- **Core single-cell workflow:** quality control, normalization, dimensionality reduction, clustering and marker detection in Seurat v5
- **Visualization:** publication-ready single-cell plots and colour palettes with SeuratExtend
- **Trajectory inference:** RNA velocity, pseudotime, cell-fate probabilities and gene trends along trajectories
- **R–Python interoperability:** running scVelo and Palantir from R through `reticulate`, converting between Seurat, loom and AnnData

## 📂 Data

- **Notebook 0:** download the [10x Genomics PBMC 3k data](https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz) and unpack it so the matrix files are in `data/pbmc3k/filtered_gene_bc_matrices/hg19/`.
- **Notebooks 1–3:** example data are downloaded automatically from [Zenodo](https://zenodo.org/records/10944066) (SeuratExtend's example datasets).

## License and credits

- Notebook 0 is adapted from the Seurat tutorial by the Satija Lab (Seurat, MIT License).
- Notebooks 1–3 are adapted from the SeuratExtend vignettes by Yichao Hua and remain under **GPL-3.0**, the SeuratExtend license.
- My own notes and code are released under the MIT License. See [`NOTICE.md`](NOTICE.md) for details.

## AI assistance

The README, repository clean-up, path fixes, the rendering of notebook 3 and the badge were prepared with the assistance of Claude Opus 5.5 (Anthropic), used through Claude Code. The study notes were written by Min-Yao Jhu. All changes were reviewed and approved by the author.

## Contact

**Min-Yao Jhu**, Senior Research Scientist, Donald Danforth Plant Science Center · [minyaojhu.com](https://www.minyaojhu.com) · [GitHub](https://github.com/MinYaoJhu)
