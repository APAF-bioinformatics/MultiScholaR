
# Verification script for LipidomicsAssayData and renamed functions
devtools::load_all(".")

# 1. Check if renamed functions exist
functions_to_check <- c(
  "generateLipidQcPlots",
  "calculateLipidMissingness",
  "getLipidInternalStandardMetrics",
  "resolveLipidDuplicateFeaturesByIntensity",
  "runLipidPerAssayRuvOptimization"
)

cat("--- Checking Renamed Functions ---\n")
for (f in functions_to_check) {
  if (exists(f)) {
    cat(sprintf("[OK] Function '%s' exists\n", f))
  } else {
    cat(sprintf("[ERROR] Function '%s' NOT found\n", f))
  }
}

# 2. Check if generic proteinMissingValueImputationLimpa is defined
cat("\n--- Checking Proteomics Generic ---\n")
if (isGeneric("proteinMissingValueImputationLimpa")) {
  cat("[OK] Generic 'proteinMissingValueImputationLimpa' exists\n")
  methods <- showMethods("proteinMissingValueImputationLimpa", printTo = FALSE)
  print(methods)
} else {
  cat("[ERROR] Generic 'proteinMissingValueImputationLimpa' NOT found\n")
}

cat("\nVerification complete.\n")
