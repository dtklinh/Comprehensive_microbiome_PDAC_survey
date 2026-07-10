---
title: "WBS — Intratumoral Microbiome Decontamination Story Line"
output:
  html_document:
    df_print: paged
---

> **Reading guide — priority tiers** [[*Essential*]]{style="color:red"} = must be in the paper for the argument to stand. [[*Better*]]{style="color:orange"} = strengthens the story, worth doing if time allows. [[*Best*]]{style="color:green"} = stretch goal / reviewer-proofing.
>
> **Critical path (Essential-only skeleton):** 1 → 1.1 → 1.2 (main) → 2.1 → 3 → 3.1 → 3.1.1. If time runs short, this subset alone still tells a complete, defensible story: NCTs were collected systematically, they carry a distinct contaminant signature that differs by control type and by lab personnel, and a longitudinal-survey-based benchmark shows decontamination reliably removes it. Everything else (Sections 1.3, 2.2, 3.1.2, 3.1.3, 3.2, 4) is corroborating or reviewer-facing evidence layered on top.

------------------------------------------------------------------------

## 1. A longitudinal contaminant survey reveals a lab-specific contamination profile. [[*Essential*]]{style="color:red"}

***Note:*** for simplicity, only Wrench normalization will be used throughout this section (state this once here so it doesn't need repeating in every subsection).

### 1.1 Negative controls (NCTs) were collected at every step of the pipeline. [[*Essential*]]{style="color:red"}

- **Figure 1A:** NCT collection scheme and associated metadata. [[*Essential*]]{style="color:red"}
- **Table 1:** NCT type × year, cells = number of NCT samples. [[*Essential*]]{style="color:red"}
- **Content:** describe the collection design and comment on the number/diversity of NCT types — this sets up why a *longitudinal* survey (rather than a single batch of controls) is necessary. [[*Essential*]]{style="color:red"}

### 1.2 Different NCT types harbor distinct contaminant profiles. [[*Essential*]]{style="color:red"}

- **Figure 1B:** alpha diversity of NCT by control type. [[*Essential*]]{style="color:red"}
- **Figure 1C:** beta diversity of NCT by control type. [[*Essential*]]{style="color:red"}
- **Figure S1:** supplementary alpha/beta diversity by year / season / lab personnel (LP). [[*Better*]]{style="color:orange"}
- **Content:**
  - Describe alpha/beta diversity by control type; highlight significant differences. [[*Essential*]]{style="color:red"}
  - Describe diversity by year/season/LP and flag notable findings — this is the bridge into Section 2. [[*Better*]]{style="color:orange"}

### 1.3 Many environmental taxa are found in NCTs, at varying abundance. [[*Better*]]{style="color:orange"}

- **Figure 1D:** heatmap (abundance & prevalence) of top genera in NCT samples.
  - Call out clusters or other visible patterns; describe a few representative taxa spanning well-known contaminants and taxa that overlap with the tumor-resident community (this overlap is the crux of why decontamination is hard — flag it explicitly, it foreshadows Section 3).
- **Figure S2:** same heatmap at other taxonomic ranks / metrics (abundance vs. prevalence).
- **Figure 1E:** differential abundance analysis (MaAsLin2, species rank) by control type.
  - Discuss insights, with appropriate caveats/speculation where the biology is uncertain.
- **Figure S3:** DAA at other taxonomic ranks and with other tools (ANCOM-BC2, ALDEx2), as robustness checks.

------------------------------------------------------------------------

## 2. Laboratory personnel influence the microbiome profile more than the processing environment does. [[*Better*]]{style="color:orange"}

*Rationale:* having shown in 1.2 that contamination varies by LP, this section asks *why* — identifying the dominant contamination source is what justifies the multi-layer statistical approach in Section 3 rather than a simpler environment-only correction.

### 2.1 LP has a stronger effect on the NCT microbiome profile than environment does. [[*Essential*]]{style="color:red"}

- **Figure 2A:** alpha diversity by LP and by environment. [[*Essential*]]{style="color:red"}
- **Figure 2B:** beta diversity, grouped by environment, colored/shaped by LP. [[*Essential*]]{style="color:red"}
- **Content:** describe and quantify the differences (effect size, not just p-values). [[*Essential*]]{style="color:red"}

### 2.2 Oral and skin commensal flora mainly drive the LP effect. [[*Best*]]{style="color:green"}

- **Figure 2C:** differential abundance by LP; describe which taxa (oral/skin genera) dominate. [[*Better*]]{style="color:orange"}
- **Figure S4:** DAA by LP/environment at other ranks and with other tools. [[*Best*]]{style="color:green"}

------------------------------------------------------------------------

## 3. A multi-layer statistical testing approach reliably removes contaminants. [[*Essential*]]{style="color:red"}

### 3.1 Benchmarking is needed, especially for low-biomass samples such as intratumoral ones. [[*Essential*]]{style="color:red"}

- Literature review and justification for benchmarking decontamination methods, with emphasis on why low-biomass samples are especially prone to false removal/retention.
- Introduce the three complementary benchmarking strategies below, used precisely because no ground truth exists.

#### 3.1.1 Benchmarking against the longitudinal negative-control survey. [[*Essential*]]{style="color:red"}

- Explain the composite contamination score conceptually, without heavy technical detail. [[*Essential*]]{style="color:red"}
- **Figure 3A:** paired Wilcoxon test, non-decontaminated (ND) vs. each decontamination method, with p-values. [[*Essential*]]{style="color:red"}
- Highlight significant differences and interpret. [[*Better*]]{style="color:orange"}
- **Supplementary table:** pairwise comparisons among ND and all decontamination methods. [[*Best*]]{style="color:green"}
- **Figure S5:** linear mixed model (random effect = SampleID) as a robustness check on the paired Wilcoxon result. [[*Best*]]{style="color:green"}

#### 3.1.2 Benchmarking with technical replicates. [[*Better*]]{style="color:orange"}

- Explain the intra-sample-distance approach using technical replicates.
- **Figure 3B:** paired Wilcoxon test with p-values, ND vs. 4 decontamination methods.
- Highlight significant differences and identify the best-performing method.
- **Supplementary table:** pairwise comparisons among ND and decontamination methods. [[*Best*]]{style="color:green"}
- **Figure S6:** LMM as an alternative/robustness check. [[*Best*]]{style="color:green"}

#### 3.1.3 Benchmarking with an overall-overlap ratio. [[*Best*]]{style="color:green"}

- Define the overlap ratio and explain its use for benchmarking.
- Describe the comparison and identify the best-performing method.
- **Figure 3C:** paired Wilcoxon test with p-values.
- **Supplementary table:** pairwise comparisons.
- **Figure S7:** LMM robustness check.

*Note: 3.1.1–3.1.3 should converge on (or explicitly reconcile disagreement about) the same "best" decontamination method — that convergence is itself a key result worth a sentence in the text.*

### 3.2 Combining technical replicates with decontamination yields a reliable microbiome profile for low-biomass samples. [[*Better*]]{style="color:orange"}

- Describe the network-based (NJ) method applied to 2 replicates, using the intersection of the two to define the reliable profile.
- **Figure 3D:** barplot of bacterial composition in each replicate and their overlap.

------------------------------------------------------------------------

## 4. FFPE samples are not reliable for intratumoral microbiome profiling. [[*Best*]]{style="color:green"}

- Report which taxa remain after decontamination in FFPE samples and compare them against the core taxa identified in Section 3.
- **[Needs expansion]** — currently a single bullet for a "Best"-tier section that closes the paper; consider adding:
  - A figure directly comparing FFPE vs. fresh-frozen (or the reference low-biomass workflow) post-decontamination taxa (overlap/Venn or barplot, consistent with Fig. 3D's style).
  - A short quantitative statement (e.g., % of core taxa retained/lost in FFPE) rather than a qualitative-only comparison.

------------------------------------------------------------------------

## 5. Synthesis / Discussion framing. [[*Better*]]{style="color:orange"} — *(new, not in original)*

Not a results section, but worth planning now so the figures above are built with the takeaway in mind:

- One paragraph tying 1–2 together: contamination is systematic and traceable to lab personnel, not random noise — this is *why* a statistical rather than ad hoc approach is required.
- One paragraph tying 3 together: triangulating three independent benchmarks (no ground truth needed) is the methodological contribution.
- One paragraph on 4: practical guidance/caveat for the field on FFPE use in microbiome studies.
- Explicit statement of the recommended decontamination method + workflow for other low-biomass microbiome studies (this is likely the single sentence reviewers will quote).
