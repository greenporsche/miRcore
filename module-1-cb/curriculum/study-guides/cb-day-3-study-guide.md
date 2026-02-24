# CB Day 3 Study Guide: Statistics and Volcano Plots

## Learning Objectives
By the end of this day, you should be able to:
- Explain the concept of a null hypothesis and why it is the starting point of statistical testing
- Interpret p-values correctly and understand what they do and do not tell you
- Describe how sample size affects statistical confidence
- Recognize a Gaussian (normal) distribution and understand standard deviation
- Calculate and interpret log2 fold change values
- Explain the purpose of the -log10 transformation for p-values
- Read and interpret a volcano plot to identify significant differentially expressed genes
- Understand the multiple testing problem and why FDR correction is necessary
- Use GEO2R to analyze your group's selected GEO dataset and export results

## Key Vocabulary
| Term | Definition |
|---|---|
| Null Hypothesis (H0) | The default assumption in a statistical test that there is no difference between groups being compared |
| Alternative Hypothesis (H1) | The claim that there IS a real difference between groups — what you are trying to find evidence for |
| P-value | The probability of observing data as extreme as (or more extreme than) what was measured, assuming the null hypothesis is true |
| Significance Threshold (alpha) | The cutoff below which p-values are considered statistically significant; conventionally set at 0.05 (5%) |
| t-test | A statistical test that compares the means of two groups (e.g., disease vs. control) to determine if the difference is statistically significant |
| Gaussian Distribution | A symmetric, bell-shaped distribution where most data points cluster around the mean; also called a normal distribution |
| Standard Deviation (SD) | A measure of how spread out data values are around the mean; approximately 68% of data falls within 1 SD of the mean |
| Sample Size (n) | The number of observations in each group; larger sample sizes produce more reliable statistical results |
| Fold Change | The ratio of a measurement in one condition to the measurement in another (Disease / Control); fold change of 2 means twice as much |
| Log2 Fold Change (log2FC) | The base-2 logarithm of the fold change; positive values = up-regulation, negative values = down-regulation, zero = no change |
| -log10(p-value) | The negative base-10 logarithm of the p-value; transforms small p-values into large positive numbers for visualization |
| Volcano Plot | A scatter plot with log2FC on the x-axis and -log10(p-value) on the y-axis; used to visualize both magnitude and significance of gene expression changes |
| Multiple Testing Problem | When performing thousands of statistical tests simultaneously (one per gene), some will produce false positive results by chance alone |
| False Discovery Rate (FDR) | A correction method for multiple testing; the FDR-adjusted p-value estimates the proportion of false positives among all results called significant |
| Adjusted P-value | A p-value that has been corrected for multiple testing (also called q-value); use this instead of raw p-value in genomics |
| Up-regulated | A gene whose expression is significantly higher in the disease group compared to the control (positive log2FC) |
| Down-regulated | A gene whose expression is significantly lower in the disease group compared to the control (negative log2FC) |
| GEO2R | An NCBI web tool that performs differential expression analysis on GEO datasets by comparing user-defined sample groups |

## Concept Summaries

### The Null Hypothesis and P-values
Statistical testing in genomics always starts with a null hypothesis: the assumption that there is no difference in gene expression between the control and disease groups. For each gene, we perform a test (typically a t-test) to see whether the observed difference is large enough to be unlikely under this assumption.

The p-value quantifies this "unlikeliness." A p-value of 0.03 means there is only a 3% probability of seeing a difference this large (or larger) if the gene's expression were actually the same in both groups. When the p-value falls below our chosen threshold (typically 0.05), we "reject" the null hypothesis and conclude the gene is differentially expressed. Crucially, the p-value is NOT the probability that the null hypothesis is true — it is the probability of the data given the null hypothesis. This distinction matters and is a common source of confusion.

### Sample Size and Confidence
Sample size profoundly affects statistical power — your ability to detect real differences. With small samples (e.g., 3 per group), even large true differences may not produce small p-values because the data is too noisy. With large samples (e.g., 100 per group), even tiny differences can become statistically significant.

This creates an important practical consideration: statistical significance does not automatically mean biological importance. A gene with a log2FC of 0.1 (only 7% expression change) might have a tiny p-value in a large study, but such a small change is probably not biologically meaningful. This is why researchers look at both p-value AND fold change together, which is exactly what a volcano plot helps you do.

### Log Transformations
Raw gene expression data and fold changes span enormous ranges, making them difficult to compare and visualize. Log transformations compress these ranges into more manageable numbers. In genomics, two log transformations are especially important.

**Log2 fold change** converts expression ratios to a symmetric scale. A gene expressed at 8x the control level has a fold change of 8 and a log2FC of 3. A gene expressed at 1/8 the control level has a fold change of 0.125 and a log2FC of -3. This symmetry is useful: a 3-fold increase and a 3-fold decrease are the same distance from zero on the log scale. A log2FC of 0 means no change at all.

**-log10(p-value)** transforms p-values for visualization. Very small p-values like 0.0001 become large positive numbers (-log10(0.0001) = 4), making them easy to see on a plot. The significance threshold of p = 0.05 corresponds to -log10(0.05) = 1.3 on this scale.

### Volcano Plots
A volcano plot is one of the most useful visualizations in genomics. It displays log2FC on the x-axis and -log10(p-value) on the y-axis, placing every tested gene as a single dot. The plot gets its name from its shape — it often resembles a volcano, with significant genes erupting upward from a base of non-significant genes.

The plot is divided into regions of interest by drawing threshold lines: vertical lines at log2FC = -1 and +1 (corresponding to 2-fold changes), and a horizontal line at the -log10(p-value) threshold for significance. Genes in the upper-right region are significantly up-regulated (high fold change AND high significance). Genes in the upper-left are significantly down-regulated. Genes in the bottom-center are neither significant nor biologically meaningful. This combined view prevents you from being misled by genes that are significant but trivially changed, or dramatically changed but not statistically reliable.

### The Multiple Testing Problem and FDR
When you test 20,000 genes for differential expression at a significance threshold of 0.05, you expect about 1,000 genes to appear "significant" purely by chance (20,000 x 0.05 = 1,000), even if no genes are truly differentially expressed. This is the multiple testing problem, and it is a fundamental challenge in genomics.

The False Discovery Rate (FDR) correction addresses this by adjusting each p-value to account for the number of tests performed. The adjusted p-value (sometimes called q-value) estimates the proportion of false positives among your significant results. When GEO2R reports an "adj.P.Val" column, this is the FDR-corrected p-value, and you should always use it instead of the raw p-value. A gene with raw p = 0.001 might have an adjusted p-value of 0.15 after correction, meaning it is not actually significant when the multiple testing context is considered.

### Using GEO2R
GEO2R is an interactive web tool provided by NCBI that allows you to analyze any GEO dataset without downloading data or writing code. You define two groups of samples (typically control and disease), assign each sample to its group, and click "Analyze." GEO2R performs statistical testing on every gene and returns a ranked table of results including gene symbols, log2FC, raw p-values, and adjusted p-values.

The results can be viewed as a "Top 250" table or exported as a full text file (TSV format) for further analysis in a spreadsheet. When exporting, the file may not open cleanly in Excel or Google Sheets — you may need to open it in a text editor first, copy the contents, and paste them into a spreadsheet. Once in a spreadsheet, you can sort and filter genes by p-value, fold change, or other criteria to generate gene lists for pathway analysis on Day 4.

## How It Connects
- **Previous Day**: On Day 2, you learned about gene expression and microarrays, and selected a GEO dataset. Today, you apply statistical methods to that dataset to separate the signal (truly differentially expressed genes) from the noise (random variation).
- **Next Day**: On Day 4, you will take your list of significant genes and investigate what they do using protein interaction networks (STRING) and biological pathway databases (KEGG). The quality of today's statistical analysis directly determines the quality of tomorrow's biological insights.

## Review Questions
1. What is the null hypothesis when comparing gene expression between a disease group and a control group?
2. Explain in plain language what a p-value of 0.001 means. What does it NOT mean?
3. Why does increasing sample size lead to more reliable (smaller) p-values for true differences?
4. Convert the following fold changes to log2 fold change: 2x, 4x, 1/2, 1/8, 1x.
5. On a volcano plot, a gene sits at coordinates (log2FC = +4, -log10 p-value = 5). Describe this gene's expression change in words.
6. If you test 10,000 genes at alpha = 0.05 and none are truly differentially expressed, how many false positives do you expect?
7. Why should you use the adjusted p-value (FDR) instead of the raw p-value when analyzing microarray results?
8. A gene has a raw p-value of 0.0005 but an adjusted p-value of 0.20. Should you consider this gene significant? Explain your reasoning.

## Further Reading
- GEO2R Tutorial: https://www.ncbi.nlm.nih.gov/geo/info/geo2r.html
- Khan Academy — P-values and Significance Tests: https://www.khanacademy.org/math/statistics-probability/significance-tests-one-sample
- StatQuest — False Discovery Rate: https://www.youtube.com/watch?v=K8LQSvtjcEo
- Nature Methods — Volcano Plots: https://www.nature.com/articles/nmeth.2837
