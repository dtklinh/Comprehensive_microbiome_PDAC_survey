## if nesscessary, install libglpk-dev, libgsl-dev, libfl-dev
## sudo apt install libcairo2-dev libfontconfig1-dev pandoc
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("phyloseq", "microbiome", "ComplexHeatmap"), update = FALSE)



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
  "VennDiagram",
  "DirichletMultinomial",
  "mia",
  "ANCOMBC",
  "emmeans",
  "pbkrtest",
  "glmnet",
  "torch",
  "devtools",
  "skimr",
  "rstatix",
  "ggpubr",
  "reshape2",
  "ALDEx2"
)

to_install <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(to_install)) BiocManager::install(to_install)

if(!"microViz" %in% rownames(installed.packages())){
  install.packages(
    "microViz",
    repos = c(davidbarnett = "https://david-barnett.r-universe.dev", getOption("repos"))
  )
}
