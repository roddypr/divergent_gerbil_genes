#!/usr/bin/R

# Performs likelihood-ratio test on PAML output,
  # comparing branch and null models

# ---------------------------------------------------------------------------- #

tests           <- read.table("results/likelihoods.txt")
colnames(tests) <- c("id", "model", "lnl", "np")

# ---------------------------------------------------------------------------- #
# Divide each model into seqparate columns
null_tests   <- tests[tests$model == "null", ]
branch_tests <- tests[tests$model == "branch", ]

stopifnot(null_tests$id == branch_tests$id)

colnames(null_tests)   <- paste(colnames(null_tests), "null", sep = "_")
colnames(branch_tests) <- paste(colnames(branch_tests), "branch", sep = "_")

tests <- cbind(null_tests, branch_tests)

# ---------------------------------------------------------------------------- #
# Test parameter
tests$D <- 2 * (tests$lnl_branch - tests$lnl_null)
# stopifnot(tests$D > 0)

# Degrees of freedom
tests$df <- tests$np_branch - tests$np_null
stopifnot(tests$df > 0)

# ---------------------------------------------------------------------------- #
# P-value
tests$p.val <- pchisq(tests$D, tests$df, lower.tail=FALSE)

# ---------------------------------------------------------------------------- #
# False Discovery Rate correction for mutiple testing
library(qvalue)
qobj <- qvalue(p = tests$p.val)

tests$q.val <- qobj$qvalues

# Plots
pdf("results/FDR.pdf")
hist(qobj)
plot(qobj)
dev.off()

# ---------------------------------------------------------------------------- #
# Parse column names to print
tests <- tests[, !colnames(tests) %in% c("id_branch", "model_null", "model_branch")]

colnames(tests)[1] <- "orthogroup"

# ---------------------------------------------------------------------------- #
# Write out
write.csv(tests, file = "results/log_ratio_test.csv", row.names=FALSE, quote=FALSE)

# ---------------------------------------------------------------------------- #
