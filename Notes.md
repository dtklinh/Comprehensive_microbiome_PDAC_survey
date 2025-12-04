# Note

1)  Normalization:

-   Wrench normalization (WN): decontam approaches should not be the condition for WN, instead using other medical caterories such as sex, treatment, condition, which will also be used for the differential analysis.
-   WN w.r.t TA or mouse gender.
-   Alternatively, applying ANCOM_BC

# Prompt

Act as biostatistician, I want to evaluate the effect of several method which are to remove contaminants. Each method produce one table, rows are 18 samples, column consists of the follows:column A is the oroginal total number of species,  column B is the number of species kept by the method (so they are significant), and column C is the number of column B overlap with given putative contaminant list. Write an R code to do the statistics to compare different methods. Combine all column into a single metric for statistical analysis

## Choose test model to distinguish true signal and contaminant
Act as bio-statistician, I have collect N1 true samples, in which n1 samples have E.coli. On the background, I randomly collect N2 sample, and in which there is n2 samples have E.coli. Design a statistical test in R, which determine if E.coli is a real species in true samples or it is just species in background