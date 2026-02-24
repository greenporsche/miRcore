# CB Day 4 Exercises: GEO2R Analysis, Pathways, and STRING

## Warm-Up Questions
1. What is the difference between looking at a single gene and looking at a pathway?
2. If you have a list of 250 differentially expressed genes, why is it more informative to find enriched biological processes than to research each gene individually?
3. What does it mean when we say a biological process is "enriched" in a gene list?

## Hands-On Exercises

### Exercise 1: Preparing Your Gene List
**Objective**: Filter and prepare a gene list from GEO2R results for pathway analysis.
**Instructions**:
1. Open your GEO2R results in Google Sheets (exported on Day 3).
2. Filter out all rows where the adjusted p-value is greater than 0.05.
3. Create three separate gene lists by copying gene symbols into new tabs:
   - **Tab 1 — Top 250 by p-value**: Sort by adjusted p-value (smallest first), copy the top 250 gene symbols
   - **Tab 2 — Up-regulated**: From the filtered list, select genes with log2FC > 0, sort by log2FC (largest first), copy the top 100 gene symbols
   - **Tab 3 — Down-regulated**: From the filtered list, select genes with log2FC < 0, sort by log2FC (most negative first), copy the top 100 gene symbols
4. Record how many total genes passed your p-value filter.

**Expected Output**: Three new tabs in your Google Sheet, each with a clean list of gene symbols (one per row, no extra columns).

### Exercise 2: STRING Network Analysis
**Objective**: Build a protein interaction network and identify key interactions.
**Instructions**:
1. Go to [STRING](https://string-db.org).
2. Click "Multiple Proteins" in the search menu.
3. Paste your top 250 gene list (Tab 1 from Exercise 1).
4. Select "Homo sapiens" as the organism and click Search.
5. On the resulting network, answer these questions:
   - How many nodes (proteins) are in the network?
   - How many edges (interactions) are shown?
   - Do you see any clusters (groups of proteins that are heavily connected to each other)?
6. Change the display to show only "Known Interactions" (cyan and magenta edges). How does the network change?
7. Save a screenshot of your network and add it to your group research document.

**Expected Output**: A STRING network screenshot with annotations noting the number of nodes, edges, and visible clusters.

### Exercise 3: GO Enrichment Analysis
**Objective**: Identify enriched biological processes in your gene list.
**Instructions**:
1. In your STRING results, click the "Analysis" tab at the bottom of the page.
2. Find the "Biological Process (Gene Ontology)" section.
3. Record the top 5 enriched biological processes:

| Rank | GO Term (Biological Process) | FDR (adjusted p-value) | Number of Genes Matching |
|---|---|---|---|
| 1 | | | |
| 2 | | | |
| 3 | | | |
| 4 | | | |
| 5 | | | |

4. For each GO term, Google the term and write a one-sentence definition in your own words.
5. For at least 2 of the GO terms, explain how they might relate to your disease of focus.
6. Export the full enrichment results: scroll to the bottom of the Analysis tab, click "Save/Export," and download the Biological Process TSV file. Import it into a new tab in your Google Sheet.

**Expected Output**: A completed table with 5 GO terms, their FDR values, student-written definitions, and a TSV imported into Google Sheets.

### Exercise 4: Comparing Up-Regulated vs. Down-Regulated Pathways
**Objective**: Discover different biological stories in up-regulated versus down-regulated genes.
**Instructions**:
1. Repeat Exercise 2 using your up-regulated gene list (Tab 2).
2. Record the top 3 enriched biological processes.
3. Repeat Exercise 2 using your down-regulated gene list (Tab 3).
4. Record the top 3 enriched biological processes.
5. Compare the two lists:
   - Are the enriched processes the same or different?
   - What biological story do the up-regulated processes tell?
   - What biological story do the down-regulated processes tell?
   - How do these stories relate to your disease?

**Expected Output**: Two separate GO enrichment tables and a written comparison paragraph.

### Exercise 5: Exploring KEGG Pathways
**Objective**: Visualize a pathway related to your disease.
**Instructions**:
1. From your STRING GO enrichment results, pick one biological process that seems highly relevant to your disease.
2. Go to [KEGG](https://www.genome.jp/kegg/).
3. Search for a related pathway (e.g., if your GO term is "apoptotic process," search for "apoptosis").
4. Find and open the pathway map.
5. Identify at least 3 genes from your top 250 list that appear in this KEGG pathway.
6. Take a screenshot of the pathway and circle or highlight the genes that are in your gene list.
7. Write a 3–4 sentence explanation of what this pathway does and why it might be relevant to your disease.

**Expected Output**: A KEGG pathway screenshot with highlighted genes and a written explanation.

## Challenge Problems

### Challenge 1: Absolute Value Analysis
**Objective**: Use absolute log2FC to find the most dramatically changed genes regardless of direction.
**Instructions**:
1. In your Google Sheet with GEO2R results, add a new column: "Absolute log2FC" = ABS(log2FC).
2. Sort by this new column (largest first).
3. Take the top 250 genes by absolute log2FC (making sure all have adj. p-value < 0.05).
4. Run this list through STRING.
5. Compare the GO enrichment results with your p-value-sorted analysis from Exercise 3.
   - What new biological processes appear?
   - What processes disappeared?
   - Which approach (sorted by p-value vs. sorted by absolute fold change) gives you more biologically interpretable results for your disease?

### Challenge 2: Network Interpretation Deep Dive
**Objective**: Go beyond the default STRING analysis to extract biological insights.
**Instructions**:
1. In your STRING network (top 250 genes), click "Clusters" to group the network.
2. Identify the 2–3 largest clusters.
3. For each cluster:
   - List the genes in the cluster
   - Run just those genes through GO enrichment (you can paste a subset into a new STRING search)
   - What is the main biological theme of each cluster?
4. Write a paragraph synthesizing your findings: How do these clusters relate to your disease? Do they represent different aspects of the disease (e.g., immune response vs. cell death vs. metabolism)?

---

## Answer Key

### Warm-Up Answers
1. A single gene provides information about one specific protein or function. A pathway shows how multiple genes and proteins work together in a coordinated biological process. Diseases typically involve entire pathways being disrupted, not just individual genes.
2. With 250 genes, researching each one individually would be overwhelming and you would miss the bigger picture. Enrichment analysis reveals the common themes and shared functions among your genes, pointing to the biological pathways most affected by the disease.
3. A biological process is "enriched" in a gene list when that process appears more frequently than you would expect by random chance. For example, if 15% of genes in the genome are involved in immune response, but 40% of your differentially expressed genes are involved in immune response, then immune response is enriched.

### Exercise Answers
1. **Gene List Preparation**: Answers vary by dataset. Key checkpoints: all genes in filtered list should have adj. p-value < 0.05; Tab 2 should contain only positive log2FC values; Tab 3 should contain only negative log2FC values.

2. **STRING Network**: Answers vary by gene list. Students should report specific numbers for nodes and edges, note whether clusters are visible, and observe that filtering to known interactions only (cyan/magenta) typically shows fewer but more reliable connections.

3. **GO Enrichment**: Answers vary by disease. Key checkpoints: FDR values should be < 0.05 for the top terms; definitions should be in the student's own words (not copy-pasted); disease connections should be specific and logical.

4. **Up vs. Down Comparison**: Answers vary. Common pattern: up-regulated genes may show enrichment for inflammatory or stress response processes, while down-regulated genes may show enrichment for normal cellular maintenance processes. The specific pattern depends on the disease.

5. **KEGG Pathways**: Answers vary. Students should select a pathway that logically connects to their GO enrichment results. Highlighted genes should actually appear in both the KEGG pathway and their gene list.

### Challenge Answers
1. **Absolute Value Analysis**: Sorting by absolute log2FC emphasizes genes with the biggest biological changes regardless of direction. This often surfaces different pathways than p-value sorting because some genes with moderate p-values have very large fold changes. Neither approach is "better" — they answer different questions. P-value sorting finds the most statistically reliable changes; fold change sorting finds the most biologically dramatic changes.

2. **Network Clusters**: Clusters represent functional modules — groups of proteins that work together. In disease research, different clusters often represent different aspects of the disease mechanism. For example, in cancer, you might see one cluster for cell cycle control, one for immune evasion, and one for metabolic reprogramming. Synthesizing cluster findings helps build a comprehensive model of disease biology.
