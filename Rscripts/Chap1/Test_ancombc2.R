## test
library(ANCOMBC)
data(dietswap, package = "microbiome")

out = ancombc2(data = pseq, tax_level = "family",
               #fix_formula = "nationality + timepoint + bmi_group",
               fix_formula = "Person + Control_type",
               rand_formula = NULL,
               p_adj_method = "holm", pseudo_sens = TRUE,
               prv_cut = 0.10, lib_cut = 500, s0_perc = 0.05,
               group = "Control_type", struc_zero = TRUE, neg_lb = TRUE,
               alpha = 0.1, n_cl = 10, verbose = TRUE,
               global = TRUE, pairwise = TRUE
               # dunnet = TRUE, 
               # trend = TRUE,
               # iter_control = list(tol = 1e-2, max_iter = 1, verbose = TRUE),
               # em_control = list(tol = 1e-5, max_iter = 1),
               # lme_control = lme4::lmerControl(),
               # mdfdr_control = list(fwer_ctrl_method = "holm", B = 1),
               # trend_control = list(contrast =
               #                        list(matrix(c(1, 0, -1, 1),
               #                                    nrow = 2,
               #                                    byrow = TRUE)),
               #                      node = list(2),
               #                      B = 1)
               )
