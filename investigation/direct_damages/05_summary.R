# =============================================================================
# Direct Damages Investigation — Summary Script
# Runs all steps in sequence.
# =============================================================================

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

message("=== Step 1: Build data ===")
source("01_data.R")

message("\n=== Step 2: Direct damage IRFs ===")
source("02_direct_irfs.R")

message("\n=== Step 3: Direct vs indirect comparison ===")
source("03_comparison.R")

message("\n=== Step 4: Heterogeneity / mediation ===")
source("04_heterogeneity.R")

message("\nAll done. Figures in figures/ — compile direct_damages_slides.tex for presentation.")
