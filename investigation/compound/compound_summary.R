# compound_summary.R
# Run all compound shock analysis scripts in order

cat("========================================\n")
cat("Compound vs Isolated Shock Analysis\n")
cat("========================================\n\n")

cat("=== Step 1: Diagnostics ===\n")
source("compound_diagnostics.R")

cat("\n=== Step 2: Estimation ===\n")
source("compound_estimation.R")

cat("\n=== Step 3: Robustness ===\n")
source("compound_robustness.R")

cat("\n=== Step 4: Figures ===\n")
source("compound_figures.R")

cat("\n========================================\n")
cat("All compound analysis complete.\n")
cat("========================================\n")
