Used R to answer the questions: 
  What happens to the fibroblast transcriptome when the cells become senescent? 
  What then happens to the fibroblast transcriptome when senescent cells have their mitochondria depleted?
  
Grade: A4

Input Required:
1. Expression Table with sample names as column names and gene IDs as row names
2. Sample Metadata
3. Annotations table
4. Differential Tables between all three conditions(with log2foldchange, p-value, p adjusted)

What the code does:
1. Parses the data, marks the significant genes above with p adjusted<0.001 and log2fold >2
2. Scaling the data for PCA, performs PCA
3. Plots PCA, MA plots and expression density plots for the whole data
4. Plots Volcano plots for each pairwise comparison
5. Compares the number of differentially expressed genes across pairwise comparisons and plots heatmaps
6. Pathway analysis for each pairwise comparison
7. Complex Analysis(All three samples compared to each other):
     Draws the foldvsfold plots for each comparison, allows for plotting of signatures from heatmaps

Note: Another file with the functions that make the plots are given in the second functions_for_analysis file.
