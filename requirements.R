## if nesscessary, install libglpk-dev
## sudo apt install libcairo2-dev libfontconfig1-dev pandoc
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("phyloseq", "microbiome", "ComplexHeatmap"), update = FALSE)

install.packages(
  "microViz",
  repos = c(davidbarnett = "https://david-barnett.r-universe.dev", getOption("repos"))
)

# requirements.R
required_pkgs <- c(
  "tidyverse",
  "vegan",
  "DESeq2",
  "Wrench",
  "ggtext",
  "ggraph",
  "DT",
  "corncob",
  "VennDiagram"
)

to_install <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(to_install)) BiocManager::install(to_install)
