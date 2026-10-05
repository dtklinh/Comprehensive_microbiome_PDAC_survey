# Reviewer 1 (Comments for the Author):

Evaluating the microbiome in low-biomass samples is notoriously challenging, and potential contamination biases have compromised several high-impact studies. To address this, the authors conducted a systematic analysis of a large cohort, incorporating negative controls and samples from multiple individuals to benchmark various decontamination protocols.

Their findings demonstrate that the NJ method yields superior decontamination performance compared to alternative methods-particularly for fresh-frozen pancreatic cancer samples, though less so for FFPE samples. Consequently, the authors recommend adopting this protocol in future studies.

The study design is straightforward, and the methodology is clearly described. The figures are clean, intuitive, and easy to follow.

Comments: In their current analysis, the authors primarily rely on species diversity and count as surrogates for decontamination quality. It would be valuable to consider whether a more species-oriented analysis could provide additional insights and definitively confirm how decontamination impacts final outputs.

To maximize the utility of this work for the broader scientific community, I strongly recommend including a clear, step-by-step pathway for implementing these decontamination protocols, even if similar guidelines have been published previously.

# Reviewer 2 (Comments for the Author):

Dang et al., assembled 337 negative controls collected over four years at three points of their wet-lab pipeline (blank paraffin scratches, lysis buffer and sterile water) and sequenced them, together with mouse tumor material, by full-length 16S rRNA gene sequencing on Oxford Nanopore. They report that the resulting contaminant profiles differ by control type, by year, by season and, most strongly, by the laboratory person who processed the sample, and that the contaminants comprise both classical environmental genera and human oral commensals. They then benchmark four decontamination strategies on fresh-frozen tumors from KPC mice, scoring them with a composite yield-and-purity metric built on the negative-control background and, orthogonally, on the agreement between technical replicates. They conclude that the Nj approach performs best on both, that it recovers a plausible gut-commensal profile from fresh-frozen tumors, and that no method recovers a comparable profile from paired FFPE material, which they interpret as showing that FFPE is unsuitable for intratumoral microbiome work in low-biomass tumors. Overall this is work that is needed, as contamination of low-biomass tissue is an important problem. While the authors have generally done a great job, there are a few issues that have to be addressed or corrected before this manuscript can be published. Above all, their entire analysis of FFPE material, and the conclusions they draw from it, should be changed (see the specific comments below).

## Major points

### 1. The FFPE analysis and all conclusions drawn from it must be revised: full-length 16S is the wrong assay for FFPE-derived DNA.

Formalin fixation fragments DNA and introduces crosslinks. A full-length (\~1,500 bp) 16S amplicon requires intact 16S molecules, which are largely depleted in FFPE-derived material. This is precisely why Nejman et al. used a 5-plex PCR of short regions spanning the 16S gene and others used relatively short amplicons. The failure to recover a bacterial profile from the FFPE samples is therefore an expected consequence of the assay, and cannot be attributed to the tissue, to contamination, or to the decontamination methods. Two further observations point the same way. FFPE material generally yields more bacterial signal than fresh-frozen tissue, not less, so the direction of the result in Figure 4A is opposite to what the field normally sees. And since some taxa are recovered better from paraffin-embedded than from fresh material, the separate clustering of FFPE from paired FF samples (lines 266-269) is expected and is not evidence that the FFPE profile is wrong or inferior. The conclusions at lines 278-283 and the summary statement at lines 386-388 ("reliably captures intratumoral microbial composition in FF samples, but not in FFPE tissues") are therefore not supported. The claim should be removed, or explicitly restricted to long-read, full-length 16S workflows. Demonstrating it properly would require a short-amplicon assay, or shotgun sequencing, on the same FFPE DNA. The current statements are wrong and should be corrected in order not to mislead the field.

### 2. Negative controls with few reads should not have been discarded.

Controls containing little or no bacterial DNA are the most informative controls in the set, and removing 75 of 337 controls biases the statistics toward more contamination. The contaminant background C, the per-control-type diversity comparisons and the variance partition in Table S3 are all therefore estimated from a contamination-enriched subset. I would suggest to repeat the key analyses with all 337 controls retained, using a model that accommodates zero and low counts rather than deleting the samples, and report how the conclusions change.

#### LD response

-   Agree to use all NCT samples

-   Raw/unadjusted alpha plot

### 3. A benchmarking study needs a ground truth, and there is none here.

There is no mock community, no dilution series, no spike-in, no positive control and no orthogonal detection of bacteria in the tumors, so no analysis in the paper can be checked against a known answer. A defined mock community (e.g. ZymoBIOMICS) carried through extraction, ideally as a dilution series into a low-biomass background, would allow per-method sensitivity and false-positive rate to be measured directly, and requires no new mice. If the authors can't add it at this point, I would at least discuss that limitation in the discussion part.

### 4. The metadata needed to judge the statistics are missing, and the design is confounded.

The number of laboratory personnel contributing to the 262-sample survey is never stated, although LP is the dominant term in the variance partition. Neither is the distribution of controls across seasons, nor the definition of "season". Without these the models in Tables S2 and S3 cannot be evaluated. Please provide the full cross-tabulation (laboratory person x year x season x control type, with counts) as a supplementary table. Table S3 is likely to show severe confounding. Table S1 already shows control types unevenly distributed across years (no PCR controls in 2021; paraffin controls 2, 5, 6, 9 and 17 across 2021-2025), and personnel presumably entered and left the laboratory over the four-year window, so LP is largely nested within year and "season" is estimated from the residual variance. The partition in Table S3 should therefore be presented as descriptive rather than causal, and the abstract (lines 44-46) tempered accordingly. The most plausible driver of the "year" effect (extraction kit and reagent lot) is not reported at all; if lot records exist they should enter the model, and if not, this belongs in the limitations.

### 5. No human PDAC sample was analyzed.

The title, the abstract (notably line 52, "Application of Nj to fresh frozen PDAC samples") and the Importance statement read as though human pancreatic tumors were profiled. All tumors are murine KPC. Please make "murine" explicit in the title and throughout the abstract.

## Minor points

### 1. Figure 1B is hard to read and does not clearly support the text.

The median is not visible for all groups, so it is difficult to tell which control type is lowest - the comparison made at lines 122-123. The group order is also counterintuitive; paraffin - buffer - PCR (or the reverse) would follow the workflow. Most importantly, the plotted data do not obviously show the highest number of observed species in the paraffin controls, as claimed at lines 128-130. Please add medians and interquartile ranges, reorder the groups, and reconcile the figure with the text.

### 2. The Shannon index result appears to contradict Figure 1B.

No significant difference in Shannon index is reported among control types (lines 131-133; Table S2, all adj. p \> 0.47), yet the buffer controls sit visibly above the others in Figure 1B. Please re-check the model and the plotted data. If the discrepancy arises because Table S2 is adjusted for LP, year and season while the figure shows raw values, state this in the legend and report the unadjusted comparison alongside.

### 3. Burkholderia appears in the buffer controls but not in the paraffin controls.

In Figure 1D, Burkholderia is essentially confined to the buffer controls. By the authors' own description (lines 125-128), paraffin controls are scraped into the same lysis buffer and carried through the same extraction and PCR, so they should contain everything the buffer controls contain, the argument used to explain why paraffin controls harbor the most taxa (lines 128-130). How do the authors explain this? Different buffer lots or aliquots for the two control types, a batch or time effect, and a compositional artifact should all be considered.

### 4. Figure 1D does not match Table S4.

The "top 15 genera" in Figure 1D omit several of the most abundant genera in Table S4 (Caldibacillus, ranked 2nd at 16.9% of reads; Fervidibacillus; Sporosarcina; Bradyrhizobium; Delftia) while including less abundant ones. Since Caldibacillus is cited as a dominant genus at lines 140-142 and is central to Figure 3D, its absence is confusing. Please state the selection criterion and reconcile with Table S4.

### 5. The laboratory-personnel experiment is narrow, and its sample numbers do not add up.

This is a small and simple experiment (two operators, two environments, water controls) carrying a full figure and a dedicated Results section, which over-weights it relative to what it can support. The numbers are also inconsistent: the text (lines 163-167) describes four sterile water samples "equally distributed", repeated in two environments, which reads as eight in total, whereas Figure 2A depicts four tubes in each of four operator-by-environment cells (sixteen) and Figure 2B appears to show roughly eight points per operator. Please state the n per operator per environment and the total, make text and figure agree, and moderate the emphasis accordingly.

### 6. It is unclear what the technical replicates replicated.

Lines 409-411 and 211-213 state that 10 FF tumors were "sequenced again as technical replicates" more than two years later by a different operator. Was DNA re-extracted from the tissue, or was the original extract re-amplified and re-sequenced? This determines whether the replicate captures extraction-derived contamination (the dominant source) or only library and sequencing variance, and therefore how much weight the inter-replicate metric can bear. Please state this in Methods.

### 7. Barcode misassignment on the ONT platform is not addressed.

Index hopping between samples on the same flow cell is a major source of apparent cross-contamination at low biomass. Were tumors and negative controls multiplexed on the same flow cells, and what is the estimated crosstalk rate (e.g. from unused barcodes)? This could by itself generate taxa shared between replicates.

### 8. Figure 2C and Table S7 disagree.

Figure 2C reports p = 0.0001 whereas Table S7 gives adj. p = 0.005 for the same comparison. Please reconcile.

### 9. The figures are supplied out of order.

In the merged file the figure cited in the text as Figure 2 appears third, after Figure 3, which makes the figure calls difficult to follow. Please check the figure order and labelling in the assembled manuscript.
