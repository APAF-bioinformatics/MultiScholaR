# Enhanced RUV-III Negative Control Optimization

## Comprehensive Guide to `findBestNegCtrlPercentage()` with Integrated K-Value Optimization

---

## Table of Contents

1. [Overview](#overview)
2. [Problem Statement](#problem-statement)
3. [Mathematical Model](#mathematical-model)
4. [Algorithm Workflow](#algorithm-workflow)
5. [Penalty Function Details](#penalty-function-details)
6. [Parameter Guidelines](#parameter-guidelines)
7. [Practical Examples](#practical-examples)
8. [Usage Instructions](#usage-instructions)
9. [Expected Outcomes](#expected-outcomes)
10. [Troubleshooting](#troubleshooting)

---

## Overview

The `findBestNegCtrlPercentage()` function automatically optimizes the selection of negative control proteins for RUV-III normalization by considering **both** separation quality in canonical correlation plots **and** the resulting k value from `findBestK()`. This prevents over-optimization towards percentages that achieve good separation but require dangerously high k values that would remove biological signal.

### Key Features

- ✅ **Automated optimization** removes subjective parameter selection
- ✅ **Composite scoring** balances separation quality vs. k value penalty
- ✅ **Integrated findBestK()** provides complete optimization in one function
- ✅ **Prevents over-correction** through mathematical k penalties
- ✅ **Preserves biological signal** while removing unwanted variation

---

## Problem Statement

### The Original Challenge

**Manual Selection Issues:**
- Subjective choice of `percentage_as_neg_ctrl`
- Time-consuming trial-and-error process
- Inconsistent results across analysts
- Difficult to justify parameter choices

**Pure Separation Optimization Issues:**
- High separation scores often correlate with high k values
- k ≥ 4 typically removes real biological signal
- Over-correction leads to loss of differential expression power
- Weird PCA patterns indicate biological signal removal

### The Solution

**Composite Optimization Approach:**
```
Composite Score = Separation Quality × (1 - K Penalty)
```

This ensures that:
- Good separation is rewarded
- High k values are penalized
- Optimal balance is achieved automatically
- Biological signal is preserved

---

## Mathematical Model

### 1. Core Components

#### A. Separation Score (Unchanged)
For each percentage tested, the separation quality is calculated using one of four metrics:

```r
# Maximum difference between All and Control correlations
separation_score = max(All_correlations - Control_correlations)
```

**Available metrics:**
- `"max_difference"`: Maximum separation across k values (default)
- `"mean_difference"`: Average separation across k values
- `"auc"`: Area under the curve of differences
- `"weighted_difference"`: K-weighted average of differences

#### B. Best K Calculation (Integrated)
For each percentage tested, the optimal k value is calculated using the existing logic:

```r
best_k = findBestK(cancorplot)
```

#### C. K Penalty Calculation (NEW)
The penalty function operates in two regimes:

**Linear Penalty Zone (k ≤ max_acceptable_k):**
```r
k_penalty = k_penalty_weight × (k - 1) / (max_acceptable_k - 1)
```

**Exponential Penalty Zone (k > max_acceptable_k):**
```r
excess_k = k - max_acceptable_k
k_penalty = k_penalty_weight + (1 - k_penalty_weight) × (1 - exp(-excess_k))
```

#### D. Composite Score
```r
composite_score = separation_score × (1 - k_penalty)
```

### 2. Default Parameters

- `k_penalty_weight = 0.5` (50% penalty weight)
- `max_acceptable_k = 3` (linear penalty up to k=3, when adaptive_k_penalty = FALSE)
- `adaptive_k_penalty = TRUE` (automatically adjusts max_acceptable_k based on sample size)

---

## Algorithm Workflow

### Step-by-Step Process

1. **Input Validation & Adaptive Setup**
   - Check data object validity
   - Validate parameter ranges
   - Calculate adaptive max_acceptable_k based on sample size (if adaptive_k_penalty = TRUE)
   - Set up progress tracking

2. **Percentage Testing Loop**
   ```r
   for each percentage in percentage_range:
       - Generate negative controls: getNegCtrlProtAnova()
       - Check minimum control count (≥5 proteins)
       - Create canonical correlation plot: ruvCancor()
       - Calculate separation score: calculateSeparationScore()
       - Calculate best k: findBestK()
       - Calculate k penalty: calculateCompositeScore()
       - Store all results
   ```

3. **Optimization**
   - Find percentage with maximum composite score
   - Extract optimal parameters and plots
   - Return comprehensive results

4. **Integration**
   - Use `best_percentage` for RUV-III normalization
   - Use `best_k` directly (no need to call `findBestK()` again)
   - Apply `best_control_genes_index` for RUV-III

---

## Penalty Function Details

### Mathematical Formulation

The penalty function is designed to:
- **Gently discourage** k values in the acceptable range (k ≤ 3)
- **Heavily penalize** k values in the dangerous range (k > 3)
- **Smoothly transition** between the two regimes

### Penalty Calculation Examples

**With default parameters (k_penalty_weight = 0.5, max_acceptable_k = 3):**

| k Value | Penalty Calculation | Penalty | Multiplier | Interpretation |
|---------|-------------------|---------|------------|----------------|
| 1 | `0.5 × (1-1)/(3-1) = 0.00` | 0.000 | 1.000 | Perfect - No penalty |
| 2 | `0.5 × (2-1)/(3-1) = 0.25` | 0.250 | 0.750 | Good - Light penalty |
| 3 | `0.5 × (3-1)/(3-1) = 0.50` | 0.500 | 0.500 | Acceptable - Moderate penalty |
| 4 | `0.5 + 0.5×(1-exp(-1)) = 0.816` | 0.816 | 0.184 | Dangerous - Heavy penalty |
| 5 | `0.5 + 0.5×(1-exp(-2)) = 0.932` | 0.932 | 0.068 | Very bad - Severe penalty |
| 6 | `0.5 + 0.5×(1-exp(-3)) = 0.975` | 0.975 | 0.025 | Terrible - Near zero score |

### Visual Representation

```
Penalty Function Shape:

Penalty
   1.0 |                    ******************
       |                ****
       |             ***
   0.5 |        ****
       |     ***
       |  ***
   0.0 |**
       +----+----+----+----+----+----+----
       1    2    3    4    5    6    7    k

       |← Linear →|←   Exponential   →|
       | Gentle   |    Aggressive     |
       | Penalty  |     Penalty       |
```

---

## Parameter Guidelines

### 1. `k_penalty_weight` (Range: 0-1)

**Controls the overall strength of k penalization:**

| Value | Behavior | Use Case |
|-------|----------|----------|
| 0.3 | Gentle penalty - favors separation quality | Large datasets with expected complex batch effects |
| 0.5 | Balanced penalty (default) | Standard proteomics experiments |
| 0.7 | Aggressive penalty - strongly favors low k | Small datasets or simple experimental designs |

### 2. `max_acceptable_k` (Range: 1-5)

**Controls the threshold for exponential penalties (when adaptive_k_penalty = FALSE):**

| Value | Behavior | Use Case |
|-------|----------|----------|
| 2 | Very conservative | Small sample sizes (n < 20), simple designs |
| 3 | Standard (default) | Typical proteomics experiments |
| 4 | Permissive | Large sample sizes (n > 50), complex designs |

### 3. `adaptive_k_penalty` (Boolean) **NEW**

**Controls whether max_acceptable_k is automatically adjusted based on sample size:**

| Value | Behavior | Use Case |
|-------|----------|----------|
| TRUE | **Automatic sample size adjustment (default & recommended)** | Most proteomics experiments - lets the function decide optimal thresholds |
| FALSE | Use fixed max_acceptable_k value | When exact reproducibility with previous results is required |

#### Adaptive K Penalty Logic

When `adaptive_k_penalty = TRUE`, the function automatically calculates the optimal `max_acceptable_k` based on your dataset's sample size:

| Sample Size | Adaptive max_k | Rationale |
|-------------|----------------|-----------|
| n < 15 | 2 | **Conservative**: Small datasets need to preserve degrees of freedom |
| n = 15-39 | 3 | **Standard**: Typical proteomics experiment range |
| n = 40-79 | 4 | **Moderate**: Larger datasets can handle more correction |
| n ≥ 80 | 5 | **Permissive**: Very large datasets have abundant statistical power |

### 4. Choosing Parameters

**Recommended Approach (adaptive penalty - default):**
```r
result <- findBestNegCtrlPercentage(
  data,
  k_penalty_weight = 0.5,
  adaptive_k_penalty = TRUE  # Default: automatically sets optimal max_acceptable_k
)
```

**Conservative Approach (for special requirements):**
```r
result <- findBestNegCtrlPercentage(
  data,
  k_penalty_weight = 0.7,          # Higher penalty
  adaptive_k_penalty = FALSE,      # Use fixed penalty
  max_acceptable_k = 2             # Very conservative threshold
)
```

**Permissive Approach (for large, complex datasets):**
```r
result <- findBestNegCtrlPercentage(
  data,
  k_penalty_weight = 0.3,          # Lower penalty
  percentage_range = seq(1, 30, by = 1),  # Wider search range
  adaptive_k_penalty = TRUE        # Let function decide max_k based on sample size
)
```

**Legacy Fixed Approach (for exact reproducibility):**
```r
result <- findBestNegCtrlPercentage(
  data,
  k_penalty_weight = 0.5,
  adaptive_k_penalty = FALSE,      # Disable adaptive behavior
  max_acceptable_k = 3             # Use exact value from previous analyses
)
```

---

## Practical Examples

### Example 1: Typical Optimization Results with Adaptive Penalty

**Input scenario:** Testing percentages 1-20% with default parameters on a medium-sized dataset (n=24)

**Sample results:**
```
Adaptive penalty enabled: Sample size = 24, Adaptive max_acceptable_k = 3
   percentage separation_score best_k composite_score num_controls valid_plot
1           4             0.040      1           0.040          178       TRUE
2           8             0.065      2           0.049          356       TRUE  ← WINNER
3          12             0.070      3           0.035          534       TRUE
4          14             0.075      5           0.005          623       TRUE
```

**Analysis:**
- Adaptive penalty set max_acceptable_k = 3 for this medium dataset
- 8% wins despite lower separation score (0.065 vs 0.075)
- k=2 is much more acceptable than k=5
- Composite score correctly identifies optimal trade-off

### Example 2: Small vs Large Dataset Adaptive Behavior

**Small Dataset (n=12) with Adaptive Penalty:**
```
Adaptive penalty enabled: Sample size = 12, Adaptive max_acceptable_k = 2
Best result: 6% (separation=0.045, k=1, composite=0.045)
Rationale: Very conservative threshold protects limited degrees of freedom
```

**Large Dataset (n=85) with Adaptive Penalty:**
```
Adaptive penalty enabled: Sample size = 85, Adaptive max_acceptable_k = 5
Best result: 15% (separation=0.080, k=4, composite=0.064)
Rationale: Permissive threshold allows higher k when statistical power supports it
```

**Same Large Dataset with Fixed Conservative Settings:**
```
Fixed penalty: max_acceptable_k = 2 (adaptive_k_penalty = FALSE)
Best result: 8% (separation=0.055, k=2, composite=0.041)
Rationale: Overly conservative for large dataset, may under-correct batch effects
```

### Example 3: Problem Detection

**Warning signs in results:**
```
   percentage separation_score best_k composite_score
1           5             0.020      1           0.020
2          10             0.025      2           0.019
3          15             0.030      3           0.015
```

**Issues indicated:**
- All composite scores are very low (< 0.05)
- May indicate poor data quality or inappropriate experimental design
- Consider different normalization approach

---

## Usage Instructions

### Basic Usage

```r
# Load your normalized protein data
# (assuming you have normalised_frozen_protein_matrix_obj)

# Run optimization with default parameters (adaptive penalty enabled by default)
optimization_result <- findBestNegCtrlPercentage(
  normalised_frozen_protein_matrix_obj,
  percentage_range = seq(1, 20, by = 1),
  separation_metric = "max_difference",
  k_penalty_weight = 0.5,
  adaptive_k_penalty = TRUE,  # Default: automatically adjusts max_acceptable_k based on sample size
  verbose = TRUE
)

# Extract optimal parameters
percentage_as_neg_ctrl <- optimization_result$best_percentage
best_k <- optimization_result$best_k
control_genes_index <- optimization_result$best_control_genes_index
cancorplot_r1 <- optimization_result$best_cancor_plot

# Update configuration
config_list$ruvParameters$percentage_as_neg_ctrl <- percentage_as_neg_ctrl

# Display results including adaptive information
cat("Optimization Results:\n")
cat("  Sample size:", optimization_result$sample_size, "samples\n")
if (!is.na(optimization_result$adaptive_k_penalty_used) && optimization_result$adaptive_k_penalty_used) {
  cat("  Adaptive max_k threshold:", optimization_result$max_acceptable_k, "\n")
} else {
  cat("  Fixed max_k threshold:", optimization_result$max_acceptable_k, "\n")
}
cat("  Best percentage:", percentage_as_neg_ctrl, "%\n")
cat("  Best k value:", best_k, "\n")
cat("  Separation score:", round(optimization_result$best_separation_score, 4), "\n")
cat("  Composite score:", round(optimization_result$best_composite_score, 4), "\n")
cat("  Number of control genes:", sum(control_genes_index, na.rm = TRUE), "\n")

# View detailed results
print(optimization_result$optimization_results)
```

### Integration with Existing Workflow

**Replace this old approach:**
```r
# OLD: Manual percentage selection
percentage_as_neg_ctrl <- 5
control_genes_index <- getNegCtrlProtAnova(data, percentage_as_neg_ctrl = 5)
cancorplot_r1 <- ruvCancor(data, ctrl = control_genes_index)
best_k <- findBestK(cancorplot_r1)
```

**With this automated approach:**
```r
# NEW: Automated optimization
optimization_result <- findBestNegCtrlPercentage(data)
percentage_as_neg_ctrl <- optimization_result$best_percentage
control_genes_index <- optimization_result$best_control_genes_index
best_k <- optimization_result$best_k  # Already calculated!
cancorplot_r1 <- optimization_result$best_cancor_plot
```

### Custom Parameter Selection

```r
# For small datasets - let adaptive penalty handle the thresholds
optimization_result <- findBestNegCtrlPercentage(
  data,
  percentage_range = seq(2, 15, by = 0.5),  # Finer resolution
  k_penalty_weight = 0.7,                   # Higher penalty
  adaptive_k_penalty = TRUE,                # Let function decide max_k (will likely be 2 for small n)
  separation_metric = "mean_difference",    # Alternative metric
  verbose = TRUE
)

# For large datasets - adaptive penalty will be more permissive
optimization_result <- findBestNegCtrlPercentage(
  data,
  percentage_range = seq(1, 25, by = 1),    # Wider range
  k_penalty_weight = 0.3,                   # Lower penalty
  adaptive_k_penalty = TRUE,                # Will set max_k = 4 or 5 for large n
  separation_metric = "auc",                # Alternative metric
  verbose = TRUE
)

# For exact reproducibility with previous analyses
optimization_result <- findBestNegCtrlPercentage(
  data,
  percentage_range = seq(1, 20, by = 1),
  k_penalty_weight = 0.5,
  adaptive_k_penalty = FALSE,               # Disable adaptive behavior
  max_acceptable_k = 3,                     # Use exact historical value
  separation_metric = "max_difference",
  verbose = TRUE
)
```

---

## Expected Outcomes

### Typical Results

**Good Optimization:**
- Best k ∈ [1, 3]
- Composite score ≥ 0.02
- Clear winner in optimization results
- Interpretable PCA plots after RUV-III

**Excellent Optimization:**
- Best k ∈ [1, 2]
- Composite score ≥ 0.05
- Multiple good options in results
- Strong separation in canonical correlation plot

### Quality Indicators

**High-Quality Results:**
```
✅ Best k ≤ 3
✅ Composite score > 0.02
✅ Separation score > 0.03
✅ > 100 control genes
✅ Clear peak in optimization curve
```

**Concerning Results:**
```
⚠️ Best k > 4
⚠️ Composite score < 0.01
⚠️ All separation scores < 0.02
⚠️ Flat optimization curve
⚠️ < 50 control genes
```

### Post-RUV-III Validation

**Check these after applying RUV-III with optimized parameters:**

1. **PCA Plot**: Should show clear biological groupings, not weird patterns
2. **RLE Plot**: Should show reduced technical variation
3. **Correlation Plot**: Should show improved within-group correlations
4. **Differential Expression**: Should detect reasonable numbers of DE proteins

---

## Troubleshooting

### Common Issues and Solutions

#### Issue 1: All Composite Scores Very Low (< 0.01)

**Possible causes:**
- Poor data quality
- Inappropriate experimental design for RUV-III
- All percentages result in high k values
- Adaptive penalty may be too conservative for your specific dataset

**Solutions:**
```r
# Try more permissive parameters with adaptive penalty
result <- findBestNegCtrlPercentage(
  data,
  k_penalty_weight = 0.3,          # Lower penalty
  adaptive_k_penalty = TRUE,       # Let function decide based on sample size
  percentage_range = seq(1, 30, by = 1)  # Wider range
)

# Or override adaptive behavior if you suspect it's too conservative
result <- findBestNegCtrlPercentage(
  data,
  k_penalty_weight = 0.3,          # Lower penalty
  adaptive_k_penalty = FALSE,      # Disable adaptive
  max_acceptable_k = 4,            # More permissive than adaptive might set
  percentage_range = seq(1, 30, by = 1)
)

# Or consider alternative normalization methods
```

#### Issue 2: Adaptive Penalty Seems Too Conservative/Permissive

**If adaptive penalty seems too conservative for your dataset:**
```r
# Override with more permissive settings
result <- findBestNegCtrlPercentage(
  data,
  adaptive_k_penalty = FALSE,      # Disable adaptive
  max_acceptable_k = 4,            # More permissive than adaptive chose
  k_penalty_weight = 0.3           # Lower penalty weight
)
```

**If adaptive penalty seems too permissive for your dataset:**
```r
# Override with more conservative settings
result <- findBestNegCtrlPercentage(
  data,
  adaptive_k_penalty = FALSE,      # Disable adaptive
  max_acceptable_k = 2,            # More conservative than adaptive chose
  k_penalty_weight = 0.7           # Higher penalty weight
)
```

#### Issue 5: Need Exact Reproducibility with Previous Results **NEW**

**Problem:** Previous analyses used fixed max_acceptable_k, but new default uses adaptive penalty

**Solution:**
```r
# Disable adaptive penalty to match historical analyses
result <- findBestNegCtrlPercentage(
  data,
  adaptive_k_penalty = FALSE,      # Critical: disable adaptive behavior
  max_acceptable_k = 3,            # Use exact value from previous analysis
  k_penalty_weight = 0.5,          # Use exact value from previous analysis
  percentage_range = seq(1, 20, by = 1)  # Use exact range from previous analysis
)
```

---

## Advanced Usage

### Custom Separation Metrics

```r
# Use AUC for datasets with complex k patterns
result <- findBestNegCtrlPercentage(
  data,
  separation_metric = "auc"
)

# Use weighted difference for emphasis on higher k values
result <- findBestNegCtrlPercentage(
  data,
  separation_metric = "weighted_difference"
)
```

### Parallel Processing (Future Enhancement)

The functional programming approach makes parallelization straightforward:

```r
# Future implementation possibility
library(furrr)
plan(multisession, workers = 4)

# Replace purrr::imap with furrr::future_imap in the source code
# for parallel percentage testing
```

### Parameter Sensitivity Analysis

```r
# Test different penalty weights
penalty_weights <- c(0.3, 0.5, 0.7)
results_list <- map(penalty_weights, ~ {
  findBestNegCtrlPercentage(
    data,
    k_penalty_weight = .x,
    verbose = FALSE
  )
})

# Compare results
map2(penalty_weights, results_list, ~ {
  cat("Penalty weight:", .x, 
      "Best %:", .y$best_percentage,
      "Best k:", .y$best_k, "\n")
})
```

---

## Conclusion

The enhanced `findBestNegCtrlPercentage()` function provides a mathematically rigorous and biologically sensible approach to RUV-III parameter optimization. The new **adaptive penalty system** automatically adjusts optimization parameters based on your dataset's sample size, ensuring optimal performance across different experimental scales while maintaining the ability to override for exact reproducibility when needed.

**Key Benefits:**
- 🎯 **Automated optimization** removes subjective decisions
- ⚖️ **Balanced scoring** considers both separation and k penalties
- 🔬 **Preserves biology** while removing technical variation
- 📊 **Improves reproducibility** across analysts and datasets
- 🚀 **Streamlines workflow** with integrated parameter selection
- 🧠 **Sample size awareness** automatically adapts to dataset characteristics
- 🔒 **Reproducibility option** maintains backward compatibility when needed

**Recommended Workflow:**
1. Use default parameters (adaptive penalty enabled) for initial optimization
2. Validate results with post-RUV-III QC plots
3. Override adaptive behavior only if results suggest it's inappropriate for your data
4. For exact reproducibility, disable adaptive penalty and use historical parameters
5. Document final parameters for reproducibility

This approach ensures that your RUV-III normalization achieves optimal technical variation removal while preserving the biological signal essential for meaningful differential expression analysis, automatically adapting to the statistical constraints and opportunities presented by your specific dataset size. 