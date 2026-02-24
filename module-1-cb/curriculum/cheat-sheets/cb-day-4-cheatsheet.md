# CB Day 4 Cheat Sheet: GEO2R Analysis, Pathways, and STRING

## Commands & Tools
| Tool/Command | What It Does | Example Usage |
|---|---|---|
| GEO2R | Identifies differentially expressed genes between sample groups | Export top 250 genes sorted by p-value |
| STRING | Protein-protein interaction network database; shows how genes/proteins relate to each other | Paste gene list into "Multiple Proteins" search |
| KEGG Pathways | Kyoto Encyclopedia of Genes and Genomes — visual maps of biological pathways | Look up a pathway mentioned in STRING results |
| Gene Ontology (GO) | Classification system for gene/protein functions organized into categories | Review Biological Process terms in STRING analysis tab |
| Google Sheets (sorting) | Sort exported gene data by different criteria | Sort by absolute value of log2FC to find most differentially expressed genes |

## Key Concepts
- **Protein-Protein Interaction (PPI) Network**: A map showing which proteins physically or functionally interact with each other
- **STRING Database**: A database of known and predicted protein interactions; displays network diagrams with colored edges indicating types of evidence
- **Known Interactions (cyan/magenta edges)**: Experimentally verified or curated interactions from databases — the most reliable type
- **Gene Ontology (GO) Terms**: Standardized vocabulary describing gene functions; divided into Biological Process, Molecular Function, and Cellular Component
- **Biological Process**: A GO category describing what a gene does in the cell (e.g., apoptosis, immune response, cell division)
- **Enrichment Analysis**: Testing whether certain functions appear more often in your gene list than expected by chance in the whole genome
- **Enrichment**: Seeing more of something than expected (e.g., more apoptosis genes in your list than in the general population)
- **Depletion**: Seeing less of something than expected
- **False Discovery Rate (FDR) in STRING**: The adjusted p-value for enrichment results; FDR < 0.05 means significant enrichment
- **Pathway Analysis**: Examining which biological pathways are represented in your gene list to understand disease mechanisms
- **KEGG Pathway**: A hand-drawn map of molecular interactions and reactions in a cell for a specific process
- **Clustering**: Grouping genes in the STRING network that are closely connected to reveal functional modules
- **Absolute Value of log2FC**: Taking |log2FC| lets you find the most differentially expressed genes regardless of direction (up or down)
- **Top 250 Genes**: A common cutoff for selecting genes to input into STRING for network analysis
- **Guest Speaker (Dr. Barmada)**: Neurologist who discussed genetic therapies, gene therapy approaches, and precision medicine applications in neurological diseases

## File Formats
- **TSV export from STRING**: Tab-separated file of GO enrichment results; download from the Analysis tab -> Save/Export
- **Network image (PNG)**: Screenshot or export of the STRING interaction network for your research document
- **Google Sheets tab**: Create separate tabs for different gene subsets (up-regulated, down-regulated, lowest p-value)

## Databases & URLs
| Resource | URL | Used For |
|---|---|---|
| STRING | https://string-db.org | Protein interaction networks and functional enrichment |
| KEGG | https://www.genome.jp/kegg/ | Biological pathway maps |
| GEO2R | https://www.ncbi.nlm.nih.gov/geo/geo2r/ | Generating gene lists for STRING input |
| Gene Ontology | http://geneontology.org | Reference for GO term definitions |

## Common Pitfalls
- Pasting genes into "Protein by Name" instead of "Multiple Proteins" in STRING — use the multiple proteins search
- Forgetting to select the correct organism in STRING (default may not be human)
- Including genes with high p-values (>0.05) in your STRING analysis — filter these out first
- Only looking at up-regulated genes — down-regulated genes also carry important information
- Not exporting the STRING network image and GO results to your shared research document
- Confusing GO Biological Process with GO Molecular Function — Biological Process describes the "what" (apoptosis), Molecular Function describes the "how" (kinase activity)
- Not Googling every GO term you find — understanding each term is essential for your final presentation

## Quick Check
1. What is the difference between enrichment and depletion?
2. Why do we input our top genes into STRING instead of looking at genes one by one?
3. In a STRING network, what do cyan and magenta edges represent?
4. If 3 out of 5 apoptosis genes in the genome appear in your top 250 differentially expressed genes, what does that suggest?
5. What should you do with every GO term listed in your STRING enrichment results?
