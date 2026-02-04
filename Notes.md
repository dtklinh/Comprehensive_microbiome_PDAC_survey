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