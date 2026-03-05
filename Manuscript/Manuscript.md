---
title: "My DPC Manuscript"
author: "Linh Dang"
bibliography: referenzen.bib
output:
  html_document:
    df_print: paged
csl: nature.csl
---

# Distangling contaminants from true intratumoural microbial signals in pancreatic ductal adenocarcinoma – A benchmarking approach

## Abstract

## Introduction

## Results

***Structure for result section***

<!--[some *blue* text]{style="color:green"}. -->

-   Contaminant species survey in the lab
    -   [Description of data collection, number of samples and taxa, significant batch effects.]{style="color:green"}
    -   [Alpha/beta analyses.]{style="color:green"}
    -   [Barplot composition / heatmap plot of bacterial composition of NCT.]{style="color:green"}
    -   Remarks
-   Hood vs. bench environments
    -   Experimetal setting
    -   Alpha/beta analyses, batch effect
    -   DAA w.r.t technicians
-   Decontamination method benchmarking
-   Intratumoral microbiome in PDAC with replicas
-   (FFPE restroration)
-   Suggested pipeline

### Contamination Survey in the Lab

#### Negative Control Samples Collection

To address the issue of contamination, we comprehensively collected a large number of negative control samples (NCT) at each stage of the wet-lab sample processing workflow (Fig. 1.1), including **93** paraffin controls, **133** buffer controls, **11** PCR controls, and **3** sequencing controls, resulting in a total of 6128 taxa. We remove sequencing control samples due to insufficient reads (≤ 750) as well as taxa with extremely low prevalence and abundance. With that criteria, we further filter out and result in **203** negative control samples, comprising **113** buffer controls, **84** paraffin controls, and **6** PCR control. That results in ~~**1775**~~**2305** taxa ~~after rarefaction as normalization~~. These negative controls were accumulated over several years, with a significant expansion in control type paraffin and buffer introduced in 2022 (Table 1 ). We also recorded several possible batch such as technicians who processed samples, year and season of processing and control types.

#### Alpha & Beta Diversity

Significant differences in microbial profiles (alpha and beta diversity) were observed across different control types (Fig. 1 & 2). Similar patterns were also found for other factors such as technician, year, and season of sequencing (in Supplementary Materials Figure 1.2, 1.3). In buffer samples, the number of species and their evenness statistically significant than the others. Notably, the association measured through Chi-suqare test among batches in sample data indicates that they are not independent (Supplementary). Thus, it is not clear if the difference in bacterial profile among control types are from their own bacterial community or due to batch effects.

#### Microbial Profile of NCT Samples

The 20 most abundant taxa in whole set of NCT represented a mixture of known environmental microbes and potential human commensals (Fig. 1.4, 1.5). For example, *Sphingomonas*, a well-known environmental taxon frequently found in hospital settings, was detected with high abundance in nearly all negative controls. In contrast, human-associated taxa such as *Veillonella parvula*, previously reported in PDAC-related studies [@McKinley2023], were present in approximately 60 percent of NCT samples but at significantly lower abundance (Fig.3). This highlights the importance of not discarding all taxa found in negative controls, but instead applying appropriate decontamination approaches [@Austin2023; @Davis2018; @Nejman2020] to systematically remove likely environmental contaminants.

Notebly, after accounting for other batch effects, we investigated the differential abundance species among control types. xxxx

**Figures and Tables**

**Table 1.** Negative control samples collected over years.

| Sample ↓ / Feature → | 2021 | 2022 | 2023 | 2024 |
|----------------------|------|------|------|------|
| Buffer               | 29   | 60   | 24   | 0    |
| Paraffin             | 5    | 57   | 11   | 12   |
| PCR                  | 0    | 3    | 3    | 0    |

![](img/Chap1/nct_exp_overview.png)

**Figure 1.1.** Above: diagram of nagative sample collection at each step of the pipeline. Below: significant batch effects asociated with negative control samples.

|  |  |  |
|----|----|----|
| ![](img/Chap1/Alpha_sampletype_obseredSpecies.png) | ![](img/Chap1/Alpha_sampletype_Shannon.png) | ![](img/Chap1/Alpha_sampletype_InvSimpson.png) |

**Figure 1.2** (A–C) Alpha diversity of NCT types w.r.t number of observed species(A), Shannon Index(B), and inversed Simpson index(C).

|  |  |
|----|----|
| ![](img/Chap1/Beta_rar_sampleType.png) | ![](img/Chap1/Beta_Wrench_SampleType.png) |

**Figure 1.3** (A-B) Beta diversity of various NCT samples, with rarefaction and Wrench normalization respectively.

|  |  |
|----|----|
| ![](img/Chap1/Heatmap_wrench_species_bysum.png) | ![](img/Chap1/Heatmap_wrench_genus_bysum.png) |
| ![](img/Chap1/Barplot_SampleType_sum.png) | ![](img/Chap1/Barplot_SampleType_prev.png) |

**Figure 1.4** Bacterial heatmap of NCT samples by species and genus (top, left to right) of top 20 most abundance taxa. Bacterial barplot of top 20 genus order by abundance and prevalence respectively (bottom, left to right).

|                                    |
|------------------------------------|
| ![](img/Chap1/NCT_composition.png) |

**Figure 1.5** Bacterial composition of nagative controls, stratified by control types.

|  |  |
|----|----|
| ![](img/Chap1/Taxatree_sampletype_season.png) | ![](img/Chap1/Taxatree_sampletype_season_key.png) |

**Figure 1.6**

### Contamination - Hood versus Bench

We were intrigued by how human commensals might enter negative control samples. To investigate this, we conducted an experiment in which negative samples were processed by two technicians under different environmental conditions, as described in Fig. XXX and detailed in the Materials and Methods section. Analysis of alpha diversity revealed no significant differences across environmental conditions, instead it shows a clear distinction among technicians (Fig. XXX). Similarly, beta diversity analysis showed that samples clustered significantly by technician (\hlred{$p < 0.001$}), but not by environmental conditions alone (\hlred{$p = 0.626$}; Fig. XXX). Importantly noted, when we accounted for technicians fluctuation, the bacterial profiles of clean and normal environments are significantly clustered, as shown in Figure XXX.

Next, we investigated the differential abundance species between clean and normal conditions. Due to the small amount number of samples as well as sparsity of NCT count tables, we applied ANCOM-BC2 [@Lin2023] and ALDEx2 [@ALDEx2] as alternatives. Regarding to the extremely low read count in NCT samples, we illustrate the DAA between different condition with relaxation of p-value threholds (ANCOM-BC2 as 0.1 and ALDEx2 as 0.2).

|                                   |                                      |
|-----------------------------------|--------------------------------------|
| ![](img/Chap2/Chap2_Alpha_TA.png) | ![](img/Chap2/Chap2_Alpha_Envir.png) |

**Figure 5.** Alpha diversity of NCT sample w.r.t environmental conditions(left) and technicians (right).

![](img/Chap2/Chap2_Beta_before_after.png)

**Figure 6.** Beta diversity cluster by conditions and technicians respectively with corresponding p-values. (Left) xxx. (Right) xxx

|  |  |
|----|----|
| ![](img/Chap2/Chap2_DAA_Envir_ANCOMBC.png) | ![](img/Chap2/Chap2_DAA_Envir_ALDEx2.png) |

**Figure 7.** Differential abundance analyses between normal (bench) and clean (hood) environments by ANCOM-BC (left) and ALDEx2 (right) respectively.

### Decontamination Methods Assessment

#### Using longitudinal survey of NCT for assessment

Due to the lack of ground truth, we employed the contaminant profile derived from the above negative sample survey as true contaminant. Here we proposed a composite score comprising two elements, yield and purity. First, yield is ...

#### Using technical replica for the assessment

We randomly select 10 PDAC fresh frozen samples from the above study, and

### Intratumoral Bacterial Profile of PDAC Fresh Frozen Samples

## Discussion

## References
