---
title: "WBS"
output:
  html_document:
    df_print: paged
---

# Work Break Down Structure for the Research Story Line

## 1. Longitudinal Contaminant Survey in the lab could reveal its unique contamination profile. [[*Essential*]]{style="color:red"}

***NOTES:*** for the sack of simplicity. I will only use Wrench normalization

### 1.1 Negative controls were collected in every step in the pipeline. [[*Essential*]]{style="color:red"}

- Figure 1.A: NCT collection and associated metadata collection. [[*Essential*]]{style="color:red"}
- Table 1: Negative control types $\times$ years, values are the number of NCT samples. [[*Essential*]]{style="color:red"}
- Content: description of the collection as well as comment on the number of NCT types. [[*Essential*]]{style="color:red"}

### 1.2 Different negative control (NCT) types harbor their own contaminant profile. [[*Essential*]]{style="color:red"}

- Figure 1.B: alpha diversity of NCT w.r.t control types. [[*Essential*]]{style="color:red"}
- Figure 1.C: beta diversity of NCT w.r.t control types. [[*Essential*]]{style="color:red"}
- Figure S1: supplementary figures of alpha/beta diversity w.r.t year/season/LP. [[*Better*]]{style="color:orange"}
- Content:
  - Describe alpha and beta diversity w.r.t control types, highlight significant differences. [[*Essential*]]{style="color:red"}
  - Similarly, describe the diversity w.r.t year/season/LP and highlight important findings. [[*Better*]]{style="color:orange"}

### 1.3 Many environmental taxa are found in NCT, with different levels of abundance. [[*Better*]]{style="color:orange"}

- Figure 1.D: heatmap (abund. & prev.) of top selected taxa (genus, by abundance) in NCT samples.
- Highlight any noticeable patterns, e.g cluster, in the plot. Also describing a few taxa, including well-known contaminant and mixed of environmental taxa and tumor-resided taxa.
- Figure S2: supplementary figures for similar heatmap plots with different taxa rank and by abund/prev.
- Figure 1.E: plots of differential abundance analyses (MaAsLin2, species rank) w.r.t control types.
- Insights of the above DAA plots, as well as certain degree of explanation and speculation.
- Figure S3: supplementary figures for DAA plots with different taxa rank as well as different softwares (ANCOM-BC2, Aldex2).

## 2. Laboratory personel influences microbiome more than processing environments. [[*Better*]]{style="color:orange"}

### 2.1 LP impacts strongly on microbiome profile of NCT. [[*Better*]]{style="color:orange"}

- Figure 2.A: alpha diversity w.r.t LP and environment. [[*Essential*]]{style="color:red"}
- Figure 2.B: beta diversity, grouped by environment and icon by LP. [[*Essential*]]{style="color:red"}
- Description the differences. [[*Essential*]]{style="color:red"}

### 2.2 Oral and skin residence flora mainly drive the impact. [[*Best*]]{style="color:green"}

- Figure 2.C: differential abundant analyses w.r.t LP, and describe. [[*Better*]]{style="color:orange"}
- Figure S4: DAA plot w.r.t LP or environment with different taxa rank and different softwares. [[*Best*]]{style="color:green"}

## 3. Multiple-layer statistical testing approach could reliably remove contaminants. [[*Essential*]]{style="color:red"}

### 3.1 Benchmarking is needed, especially when working with low-biomass sample such as intratumoral one. [[*Essential*]]{style="color:red"}

- Literature review and justification for the need of decontamination method benchmarking. Emphasize on the tricky situation dealing with low-biomass samples.
- Highlight the following three different ways for benchmarking when lacking of the ground truth.

#### 3.1.1 Benchmarking with longitudinal negative survey. [[*Essential*]]{style="color:red"}

- Explanation of the composite score without much technical details. [[*Essential*]]{style="color:red"}
- Figure 3.A: paired Wilcox-test plot with p-values. [[*Essential*]]{style="color:red"}
- Highlight the significant differences and describe the results. [[*Better*]]{style="color:orange"}
- Table for supplementary: pairwise comparison among non-decontaminated (ND) data and data from decontamination methods. [[*Best*]]{style="color:green"}
- Figure S5: plot in supplementary using Linear Mixed Model with random effect as SampleID as an alternative to pair wilcox-test. [[*Best*]]{style="color:green"}

#### 3.1.2 Benchmarking with technical replicas. [[*Better*]]{style="color:orange"}

- Explanation of using the technical replicas for benchmarking, intra-sample distance.
- Figure 3.B: paired Wilco-test plot with p-values, among NP and 4 decontamination methods.
- Highlight the significant differences and emphasize on the best perfomced method.
- Table for supplementary: similarly, pairwise comparison among non-decontaminated (ND) data and data from decontamination methods. [[*Best*]]{style="color:green"}.
- Figure S6: plot in supplementary using LMM as an alternative. [[*Best*]]{style="color:green"}

#### 3.1.3 Benchmarking with overall-overlap ratio. [[*Best*]]{style="color:green"}

- Explanation what is that ratio and how it could be used for the benchmarking.
- Describe the comparison and highlight best method.
- Figure 3.C: paired Wilcox-test with p-values.
- Pairwise comparison table for supplementary.
- Figure S7 for supplementary, using LMM.

### 3.2 Technical replicas combined with deconmatimantion method could be useful to reveal reliable microbiome profile for low-biomass samples. [[*Better*]]{style="color:orange"}

- Describe using Nj-based method on 2 replica and the intersection of them to reveal the reliable microbiome profile.
- Figure 3.D: a barplot showing bacterial composition of 2 replicas and their overlap.

## 4. FFPE samples are not reliable for intratumoral microbiome profile. [[*Best*]]{style="color:green"}

- Report the remaining taxa after decontaminated and compared them with core taxa found from the above section.
