# CB Day 4 Study Guide: GEO2R Analysis, Pathways, and STRING

## Learning Objectives
By the end of this day, you should be able to:
- Export and organize GEO2R results into a spreadsheet for analysis
- Filter gene lists by adjusted p-value and sort by different criteria (p-value, log2FC, absolute log2FC)
- Use the STRING database to build protein-protein interaction networks
- Interpret STRING network diagrams, including edge colors and clusters
- Perform Gene Ontology (GO) enrichment analysis using STRING
- Explain what enrichment and depletion mean in the context of gene list analysis
- Export and interpret GO biological process results
- Navigate KEGG pathway maps and identify relevant genes
- Describe how guest speaker Dr. Barmada connected genetic therapies to precision medicine

## Key Vocabulary
| Term | Definition |
|---|---|
| STRING | Search Tool for the Retrieval of Interacting Genes/Proteins — a database of known and predicted protein-protein interactions |
| Protein-Protein Interaction (PPI) | A physical or functional association between two or more proteins; proteins rarely work alone |
| Network Node | A circle in a STRING network representing a single protein/gene |
| Network Edge | A line connecting two nodes in a STRING network, representing an interaction; different colors indicate different types of evidence |
| Known Interaction (cyan edge) | An interaction verified by curated databases |
| Known Interaction (magenta edge) | An interaction determined by laboratory experiments |
| Gene Ontology (GO) | A standardized system for classifying gene and protein functions into three categories: Biological Process, Molecular Function, and Cellular Component |
| Biological Process | A GO category describing what the gene does at the cellular level (e.g., apoptosis, immune response, DNA repair) |
| Molecular Function | A GO category describing the biochemical activity of the gene product (e.g., kinase activity, DNA binding) |
| Cellular Component | A GO category describing where in the cell the gene product is located (e.g., nucleus, membrane, cytoplasm) |
| Enrichment | The observation that a particular function or category appears more frequently in a gene list than expected by random chance |
| Depletion | The observation that a particular function or category appears less frequently than expected |
| False Discovery Rate (FDR) in enrichment | The adjusted p-value for enrichment results; FDR < 0.05 indicates statistically significant enrichment |
| KEGG (Kyoto Encyclopedia of Genes and Genomes) | A database of biological pathway maps showing how molecules interact and react within cells |
| Pathway | A series of molecular interactions and reactions in a cell that lead to a particular biological outcome |
| Clustering | Grouping proteins in a network based on their connectivity; clusters often represent functional modules |
| TSV (Tab-Separated Values) | A text file format where columns are separated by tabs; used for STRING and GEO2R exports |
| Gene Therapy | A therapeutic approach that modifies or replaces faulty genes to treat or prevent disease |

## Concept Summaries

### From Gene List to Biological Meaning
After Day 3, you have a spreadsheet full of differentially expressed genes sorted by p-value and fold change. But a list of gene names, no matter how long, does not by itself tell you what is happening biologically. The crucial next step is to ask: what do these genes have in common? Are they involved in the same biological processes? Do they interact with each other?

This is where pathway and network analysis come in. Instead of researching each of your 250+ genes individually (which would be overwhelming and miss the bigger picture), you use computational tools to find patterns. STRING identifies which of your genes encode proteins that interact with each other, and Gene Ontology enrichment testing reveals which biological functions are over-represented in your gene list. Together, these analyses transform a list of gene names into a story about which biological processes are disrupted in your disease.

### STRING: Protein Interaction Networks
STRING (string-db.org) is one of the most widely used tools in bioinformatics for understanding relationships between genes. When you input a list of genes (using the "Multiple Proteins" search), STRING creates a network diagram where each gene is a node (circle) and connections between genes are edges (lines). The colors of edges indicate the type of evidence supporting the interaction: cyan and magenta represent known interactions from curated databases and experiments, which are the most reliable.

When interpreting a STRING network, look for clusters — dense groups of interconnected proteins. These clusters often represent functional modules: groups of proteins that work together in a specific biological process. For example, a cluster of immune-related proteins in your disease network might suggest that the immune system is heavily involved in your disease. The density and structure of the network itself is informative: a highly connected network suggests your differentially expressed genes are functionally related, not random.

### Gene Ontology and Enrichment Analysis
Gene Ontology (GO) provides a standardized vocabulary for describing what genes do. When STRING performs enrichment analysis, it asks: among the ~20,000 human genes, what proportion are involved in, say, "apoptotic process"? And among YOUR top 250 differentially expressed genes, what proportion are involved in apoptotic process? If the proportion in your gene list is significantly higher than in the genome overall, that process is said to be "enriched."

Think of it like this: if you randomly selected 250 people from the general population, you would expect a certain proportion to be left-handed (about 10%). But if you selected 250 people from a left-handed bowling league, the proportion of left-handers would be much higher — left-handedness is "enriched" in your selected group. Similarly, when a biological process is enriched in your differentially expressed genes, it suggests that process is specifically affected by the disease, not a coincidence. The FDR value tells you how confident you can be in this enrichment.

### KEGG Pathways: Visualizing Biological Processes
While GO terms tell you WHAT functions are enriched, KEGG pathway maps show you HOW the molecules work together. KEGG provides hand-curated diagrams of biological pathways, showing each protein as a box and the interactions between them as arrows. These maps look like circuit diagrams for the cell.

When you look up a KEGG pathway related to one of your enriched GO terms, you can visually identify where your differentially expressed genes fit within the broader biological process. For instance, if "apoptosis" is enriched in your results, you can open the KEGG apoptosis pathway and find which specific steps in the apoptotic process are affected by the genes in your list. This provides a mechanistic understanding of the disease that goes far beyond simply knowing which genes are up or down.

### Guest Speaker: Dr. Barmada on Genetic Therapies
Dr. Sami Barmada, a neurologist specializing in memory disorders, spoke about how genetic research translates into real therapies. His work on neurodegenerative diseases illustrates the pipeline from gene discovery to treatment: first identifying the genes involved (as students did with OMIM and GEO2R), then understanding the pathways affected (as students did with STRING and KEGG), and finally developing targeted interventions.

Dr. Barmada discussed gene therapy — a cutting-edge approach that modifies faulty genes directly, either by correcting mutations, silencing harmful genes, or adding functional copies of missing genes. He highlighted that precision medicine is not just theoretical — real patients are receiving genetically targeted treatments today. This perspective reinforces why the bioinformatics skills students learn during camp have direct real-world applications in medicine and drug development.

### Practical Tips for STRING Analysis
When using STRING, keep several practical points in mind. First, always use "Multiple Proteins" search (not "Protein by Name," which is for single-protein lookups). Second, make sure to select the correct organism — the default may not be human. Third, filter your gene list before inputting it: only include genes with adjusted p-value < 0.05 to avoid introducing noise. Fourth, try analyzing different subsets: all significant genes, only up-regulated, only down-regulated, and sorted by absolute log2FC. Different subsets may reveal different biological stories. Finally, always export your results — both the network image and the GO enrichment TSV file — to your shared research document.

## How It Connects
- **Previous Day**: On Day 3, you learned statistics, performed GEO2R analysis, and generated a list of differentially expressed genes. Today, you take that gene list and discover its biological meaning through network and pathway analysis.
- **Next Day**: On Day 5, you will synthesize everything from this week into a research presentation. The STRING networks, GO enrichment results, and KEGG pathways from today will be key figures in your slides. You will also discuss the social impact of your disease, connecting the molecular findings to real-world implications.

## Review Questions
1. Why is it more informative to analyze a gene list as a network rather than studying each gene individually?
2. In a STRING network, what do cyan and magenta edges represent, and why are they considered the most reliable?
3. Explain the concept of enrichment using a real-world analogy (not the one from class).
4. You run STRING on your top 250 genes and find that "inflammatory response" has an FDR of 0.0003 and "signal transduction" has an FDR of 0.15. Which is significantly enriched? What does the FDR value tell you?
5. What is the difference between a GO Biological Process and a KEGG Pathway?
6. Describe the steps to export GO enrichment results from STRING into Google Sheets.
7. Why might analyzing only up-regulated genes give different enrichment results than analyzing only down-regulated genes?
8. How does the research workflow you have followed this week (OMIM -> GEO -> GEO2R -> STRING -> KEGG) mirror what professional bioinformaticians do?

## Further Reading
- STRING Database: https://string-db.org
- STRING Help and Documentation: https://string-db.org/cgi/help
- KEGG Pathway Database: https://www.genome.jp/kegg/pathway.html
- Gene Ontology Resource: http://geneontology.org
- Nature Protocols — STRING Analysis: https://www.nature.com/articles/nprot.2016.018
