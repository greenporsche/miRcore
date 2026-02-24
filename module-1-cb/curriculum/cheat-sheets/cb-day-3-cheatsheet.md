# CB Day 3 Cheat Sheet: Statistics and Volcano Plots

## Commands & Tools
| Tool/Command | What It Does | Example Usage |
|---|---|---|
| GEO2R | NCBI tool that compares groups of samples in a GEO dataset to find differentially expressed genes | Select control vs. disease samples, click "Analyze" |
| Interactive P-value Simulator | Visual tool to explore how sample size and effect size influence p-values | Adjust sliders to see how distributions change |
| Google Sheets (sorting/filtering) | Sort gene lists by p-value or log fold change columns | Right-click column header -> Sort Sheet A to Z |

## Key Concepts
- **Null Hypothesis (H0)**: The default assumption that there is no difference between groups (e.g., no difference in gene expression between disease and control)
- **P-value**: The probability of observing your results if the null hypothesis were true; lower p-value = stronger evidence against H0
- **Significance Threshold**: Typically p < 0.05 (5% chance the result is due to random chance); in genomics, often stricter (p < 0.01 or lower)
- **Confidence Level**: 1 minus the p-value; a p-value of 0.05 corresponds to 95% confidence
- **t-test**: A statistical test that compares the means of two groups to determine if they are significantly different
- **Gaussian (Normal) Distribution**: The bell-shaped curve that describes how data points distribute around a mean
- **Standard Deviation**: A measure of how spread out data points are from the mean; larger SD = more variability
- **Sample Size**: The number of observations in each group; larger samples give more reliable results and smaller p-values
- **Log Transformation**: Converting values to logarithmic scale to make large ranges of data easier to compare and visualize
- **Log2 Fold Change (log2FC)**: The log base 2 of the ratio of expression between disease and control; positive = up-regulated, negative = down-regulated
- **Fold Change**: How many times more (or less) a gene is expressed in disease vs. control (e.g., 2-fold = twice as much)
- **-log10(p-value)**: Negative log base 10 of the p-value; transforms small p-values into large positive numbers for plotting
- **Volcano Plot**: A scatter plot with log2FC on the x-axis and -log10(p-value) on the y-axis; helps identify genes that are both statistically significant AND biologically meaningful
- **Multiple Testing Problem**: When testing thousands of genes, some will appear significant by chance; requires correction
- **False Discovery Rate (FDR)**: An adjusted p-value that accounts for multiple testing; also called adjusted p-value or q-value
- **Up-regulated Gene**: A gene expressed at higher levels in the disease group compared to control (positive log2FC)
- **Down-regulated Gene**: A gene expressed at lower levels in the disease group compared to control (negative log2FC)

## File Formats
- **GEO2R output table**: Tab-delimited results with columns for gene ID, log2FC, p-value, adjusted p-value
- **TSV (Tab-Separated Values)**: GEO2R exports data in this format; open with Notepad (Windows) or TextEdit (Mac), then paste into Google Sheets

## Databases & URLs
| Resource | URL | Used For |
|---|---|---|
| GEO2R | https://www.ncbi.nlm.nih.gov/geo/geo2r/ | Comparing sample groups within a GEO dataset to find differentially expressed genes |
| NCBI GEO | https://www.ncbi.nlm.nih.gov/geo/ | Finding the dataset to analyze with GEO2R |

## Common Pitfalls
- Confusing p-value with probability of the hypothesis being true — p-value is the probability of the DATA given H0, not vice versa
- Forgetting that log2FC of 0 means no change, not zero expression
- Ignoring the multiple testing problem — a raw p-value of 0.05 is not reliable when testing 20,000 genes simultaneously
- Misreading volcano plots — the most interesting genes are in the upper-left (down-regulated, significant) and upper-right (up-regulated, significant) corners
- Trouble opening TSV files — if your computer opens the GEO2R export as garbled text, open it with a plain text editor first, then copy-paste into Google Sheets
- Forgetting to assign samples to the correct groups (control vs. disease) in GEO2R before running the analysis

## Quick Check
1. What does a p-value of 0.01 mean in plain language?
2. If a gene has a log2 fold change of -3, is it up-regulated or down-regulated, and by approximately how much?
3. Why do we use -log10(p-value) instead of raw p-value on volcano plots?
4. What is the False Discovery Rate and why is it needed in genomics?
5. On a volcano plot, where would you find the most important genes to study?
