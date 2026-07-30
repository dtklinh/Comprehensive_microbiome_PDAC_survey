# Note

1)  Normalization:

-   Wrench normalization (WN): decontam approaches should not be the condition for WN, instead using other medical caterories such as sex, treatment, condition, which will also be used for the differential analysis.
-   WN w.r.t TA or mouse gender.
-   Alternatively, applying ANCOM_BC

# Prompt

Act as biostatistician, I want to evaluate the effect of several method which are to remove contaminants. Each method produce one table, rows are 18 samples, column consists of the follows:column A is the oroginal total number of species,  column B is the number of species kept by the method (so they are significant), and column C is the number of column B overlap with given putative contaminant list. Write an R code to do the statistics to compare different methods. Combine all column into a single metric for statistical analysis

## Choose test model to distinguish true signal and contaminant
Act as bio-statistician, I have collect N1 true samples, in which n1 samples have E.coli. On the background, I randomly collect N2 sample, and in which there is n2 samples have E.coli. Design a statistical test in R, which determine if E.coli is a real species in true samples or it is just species in background

## Promt to use permutation test as an alternative for linear mixed model 
act as Biostatistician and give me an toy examplein in R. I have one set of samples, undergoing on different treatments, let say A and B. Using linear mixed model, I could consider that treatment as   fixed effect (I want to compare) and sample_id as random effect (Accounting for biological variation between samples). Now, I want to use the permutation test. Design the alternative permutation test which has similar performance as above linear mixed model.

## Promt to work with replica

I have a microbiome dataset including 10 tumor samples. Each sample can be sequenced many times, we call that technical replica. We assume that replicas from the sample sample should share similar bacterial profile when contaminant are removed. 
I have 3 different methods, named M1, M2, M3. They are methods to remove contamination in the original dataset. Now, after one year, from the original sample set, I acquire a set of technical replica of the original set. I ran 3 decontamination methods on those replica as well. 
Act as bistatistician, design a framework to evaluate 3 above methods M1, M2, M3. Assess them of how good that they remove contaminants. 

## Paired wilcoxon vs linear mixed effects model with random intercept - the 2 interpretations
https://stats.stackexchange.com/questions/522264/paired-wilcoxon-vs-linear-mixed-effects-model-with-random-intercept-the-2-inte


here is the SCRuB paper: https://doi.org/10.1038/s41587-023-01696-w
## Note on SCRuB
- its mechanism: reduce abundance of the contaminant taxa by subtract reads from the pooled background, but not totally eliminate the taxa. 
- SCRuB best performce with well-to-well spatial information.
-assessment metric: mainly based on present/absent information, thus read abundance adjustment in SCRuB may not be suitable for the metric. 
- However, the overl bacterial composition is comparable to Nj, suggesting its effectiveness.

## Goal:
write a short paragraph for discussion of speculation that SCRuB may not performance well.

-->
Despite SCRuB's strong overall performance in benchmark validations where the method outperformed state-of-the-art decontamination approaches by an average of 15–20 times across multiple ecosystems, data types, and sequencing depths, its relatively modest performance in our present/absent-based assessment may reflect a mismatch between the algorithm's underlying mechanism and the metric used to evaluate it. Rather than eliminating contaminant taxa outright, SCRuB operates by probabilistically subtracting reads attributable to a pooled background contamination source, leaving residual low-level abundances of contaminant taxa rather than removing them entirely. A presence/absence-based metric would therefore fail to credit this partial, abundance-based correction, since a taxon that has been substantially down-weighted but not fully zeroed out would still register as "present" and be penalized as if no decontamination had occurred. Furthermore, SCRuB is reported to achieve its best performance when well-to-well spatial (plate-position) information is available, allowing it to model cross-contamination via leakage between neighboring wells as it incorporates information regarding the spatial position of samples during processing to detect cross-contamination; the absence or incompleteness of such spatial metadata in our cohort may have limited SCRuB from reaching its full decontamination potential. Nonetheless, the overall bacterial community composition recovered by SCRuB remained comparable to that of Nj, suggesting that despite its apparent underperformance on our binary metric, the method still captured biologically meaningful signal and preserved the broader compositional structure of the samples. 

-->
SCRuB's relatively modest performance under our present/absent-based assessment may stem from a mismatch between its mechanism and the metric itself. 
Rather than eliminating contaminant taxa outright, SCRuB probabilistically subtracts reads attributable to a pooled background source, leaving residual abundances of contaminant taxa rather than removing them entirely — a partial correction that a binary presence/absence metric cannot credit, since a down-weighted but non-zero taxon still registers as "present." Performance may be further limited by the absence or incompleteness of well-to-well spatial (plate-position) metadata in our cohort, as SCRuB is reported to achieve its best results when such information is available to model cross-contamination via leakage between neighboring wells. Nonetheless, the overall bacterial composition recovered by SCRuB remained comparable to that of Nj, suggesting that despite its apparent underperformance on our metric, the method still captured biologically meaningful signal.