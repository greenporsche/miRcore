# CB Day 3 Exercises: Statistics and Volcano Plots

## Warm-Up Questions
1. What is a null hypothesis, and why do we start with one in a scientific study?
2. If you flip a coin 10 times and get 7 heads, is that surprising? What if you flip it 10,000 times and get 7,000 heads?
3. What does it mean for a result to be "statistically significant"?

## Hands-On Exercises

### Exercise 1: Understanding P-values
**Objective**: Build intuition for what p-values mean.
**Instructions**:
For each scenario below, estimate whether the p-value would be very small (< 0.01), moderate (0.01–0.05), or large (> 0.05):

| Scenario | Your P-value Estimate |
|---|---|
| A. You compare height between 5 basketball players and 5 random people. The basketball players are taller on average. | |
| B. You compare height between 500 basketball players and 500 random people. The basketball players are taller on average. | |
| C. You compare test scores between two classes of 30 students each. The averages are 82.1 and 82.3. | |
| D. You compare gene expression of TP53 between 10 cancer samples and 10 healthy samples. Cancer samples show 8x lower expression. | |

For each, explain your reasoning in one sentence.

**Expected Output**: Reasonable estimates with explanations referencing sample size and effect size.

### Exercise 2: Log Fold Change Calculations
**Objective**: Practice converting between fold change and log2 fold change.
**Instructions**:
Complete the following table:

| Gene | Control Expression | Disease Expression | Fold Change (Disease/Control) | log2(Fold Change) | Up or Down Regulated? |
|---|---|---|---|---|---|
| Gene A | 100 | 400 | | | |
| Gene B | 800 | 100 | | | |
| Gene C | 250 | 250 | | | |
| Gene D | 50 | 1600 | | | |
| Gene E | 1200 | 300 | | | |

Hint: log2(1) = 0, log2(2) = 1, log2(4) = 2, log2(8) = 3, log2(0.5) = -1, log2(0.25) = -2, log2(0.125) = -3

**Expected Output**: Completed table with correct fold changes, log2 values, and regulation direction.

### Exercise 3: Reading a Volcano Plot
**Objective**: Interpret a volcano plot to identify significant genes.
**Instructions**:
Consider a volcano plot where the x-axis is log2(Fold Change) and the y-axis is -log10(p-value). The significance thresholds are |log2FC| > 1 and -log10(p-value) > 1.3 (which corresponds to p < 0.05).

The following genes are plotted at these coordinates:

| Gene | log2FC | -log10(p-value) | Quadrant |
|---|---|---|---|
| Gene 1 | +3.5 | 4.0 | |
| Gene 2 | -2.0 | 0.8 | |
| Gene 3 | +0.3 | 5.0 | |
| Gene 4 | -4.0 | 6.0 | |
| Gene 5 | +0.1 | 0.2 | |

1. For each gene, determine: Is it significantly differentially expressed? Why or why not?
2. Which genes would you prioritize for further research?
3. What does Gene 3 represent — is it useful or misleading?

**Expected Output**: Gene 1 and Gene 4 are significantly differentially expressed (pass both thresholds). Gene 2 has large fold change but is not statistically significant. Gene 3 is statistically significant but has very small fold change (biologically not meaningful). Gene 5 is neither significant nor biologically meaningful.

### Exercise 4: Running GEO2R
**Objective**: Perform a differential expression analysis using GEO2R.
**Instructions**:
1. Go to your group's selected GEO dataset page (from Day 2).
2. Click "Analyze with GEO2R."
3. Define two groups:
   - Group 1: Control samples (label them "Control")
   - Group 2: Disease samples (label them "Disease")
4. Assign each sample to the correct group.
5. Click "Top 250" to see the results.
6. Record the top 5 genes with their:
   - Gene symbol
   - log2 fold change
   - Adjusted p-value
   - Whether they are up- or down-regulated
7. Export the full results table and paste it into a Google Sheet.

**Expected Output**: A Google Sheet tab with the GEO2R results and a written summary of the top 5 genes.

### Exercise 5: Volcano Plot Sketching
**Objective**: Create a hand-drawn volcano plot from your GEO2R data.
**Instructions**:
Using your GEO2R export:
1. Draw axes: x-axis = log2FC, y-axis = -log10(adj. p-value).
2. Draw dashed lines for your significance thresholds: vertical lines at log2FC = -1 and +1, horizontal line at -log10(0.05) = 1.3.
3. Plot at least 10 genes from your results (use dots for each gene).
4. Label any genes that fall in the "significant and large fold change" regions (upper left and upper right).
5. Color code: red for up-regulated significant, blue for down-regulated significant, gray for non-significant.

**Expected Output**: A hand-drawn or digitally sketched volcano plot with labeled significant genes.

## Challenge Problems

### Challenge 1: The Multiple Testing Problem
**Objective**: Understand why raw p-values are unreliable when testing many genes.
**Instructions**:
1. Suppose you test 20,000 genes for differential expression using a p-value threshold of 0.05.
2. If NONE of the genes were truly differentially expressed (all null hypotheses are true), how many genes would you expect to appear "significant" by chance alone? Show your calculation.
3. This is why we use the adjusted p-value (FDR). If the FDR-adjusted p-value for a gene is 0.03, what does that mean in plain language?
4. Your GEO2R results include a column called "adj.P.Val." A gene has raw p-value = 0.001 but adj.P.Val = 0.12. Should you call this gene significant? Why or why not?

### Challenge 2: Effect Size vs. Statistical Significance
**Objective**: Distinguish between statistical significance and biological importance.
**Instructions**:
Consider two genes:
- **Gene X**: log2FC = 0.15, adjusted p-value = 0.0001 (from a study with 500 samples per group)
- **Gene Y**: log2FC = 3.2, adjusted p-value = 0.08 (from a study with 4 samples per group)

1. Which gene is statistically significant? Which is biologically meaningful?
2. Why might Gene X have such a small p-value despite a tiny fold change?
3. Why might Gene Y fail the significance threshold despite a huge fold change?
4. If you had to choose one gene for follow-up research, which would you pick and why?

---

## Answer Key

### Warm-Up Answers
1. The null hypothesis is the default assumption that there is no difference between groups. We start with it because science works by trying to disprove the null — if we can show the data is very unlikely under H0, we have evidence that a real difference exists.
2. 7/10 heads (70%) is not very surprising for a small sample — random variation is high. But 7,000/10,000 heads (70%) is extremely surprising because with a large sample, random variation is much smaller. This demonstrates how sample size affects our confidence.
3. A result is statistically significant when the p-value is below a chosen threshold (usually 0.05), meaning there is less than a 5% probability that the observed difference occurred by chance.

### Exercise Answers
1. **P-value Estimates**:
   - A: Large (> 0.05) — small sample size makes it hard to be confident even though the difference is real
   - B: Very small (< 0.01) — large sample size with a real, consistent difference
   - C: Large (> 0.05) — the difference between 82.1 and 82.3 is tiny and likely due to chance
   - D: Very small (< 0.01) — 8x difference is a huge effect size and 10 samples per group provides reasonable power

2. **Log Fold Change Calculations**:

| Gene | Control | Disease | Fold Change | log2FC | Direction |
|---|---|---|---|---|---|
| Gene A | 100 | 400 | 4.0 | +2.0 | Up-regulated |
| Gene B | 800 | 100 | 0.125 | -3.0 | Down-regulated |
| Gene C | 250 | 250 | 1.0 | 0.0 | No change |
| Gene D | 50 | 1600 | 32.0 | +5.0 | Up-regulated |
| Gene E | 1200 | 300 | 0.25 | -2.0 | Down-regulated |

3. **Volcano Plot Interpretation**:
   - Gene 1 (log2FC=+3.5, -log10p=4.0): Significantly up-regulated — passes both thresholds
   - Gene 2 (log2FC=-2.0, -log10p=0.8): Large fold change but NOT significant (p > 0.05) — could be noise
   - Gene 3 (log2FC=+0.3, -log10p=5.0): Statistically significant but very small fold change — biologically not meaningful, likely an artifact of large sample size
   - Gene 4 (log2FC=-4.0, -log10p=6.0): Significantly down-regulated — passes both thresholds, strongest candidate
   - Gene 5 (log2FC=+0.1, -log10p=0.2): Neither significant nor meaningful — background noise
   - Priority genes: Gene 4 (strongest signal) and Gene 1

4. **GEO2R Analysis**: Answers will vary by group and dataset.

5. **Volcano Plot Sketching**: Evaluated on correct axis labels, threshold lines, gene placement, and color coding.

### Challenge Answers
1. **Multiple Testing**: 20,000 genes x 0.05 = 1,000 genes expected to appear significant by chance alone. An FDR of 0.03 means that among all genes called significant, approximately 3% are expected to be false positives. For the gene with raw p = 0.001 but adj.P.Val = 0.12, you should NOT call it significant because after correcting for multiple testing, the adjusted p-value exceeds 0.05.

2. **Effect Size vs. Significance**: Gene X is statistically significant (adj. p = 0.0001) but biologically trivial (0.15 log2FC = only 11% expression change). Gene Y is biologically meaningful (3.2 log2FC = ~9x change) but not statistically significant (p = 0.08). Gene X's tiny p-value comes from the massive sample size (500 per group) — even tiny differences become "significant" with enough data. Gene Y fails significance because 4 samples per group gives very low statistical power. For follow-up, Gene Y is the better choice — it shows a dramatic biological effect and likely just needs more samples to confirm. A biologically trivial difference, no matter how statistically significant, is rarely worth pursuing.
