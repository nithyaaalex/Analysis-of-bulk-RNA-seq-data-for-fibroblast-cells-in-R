#### FUNCTIONS ####
#parses the de file, annotated it, adds some columns and returns it
parse_de = function(de, em_annotated, p_threshold, fold_threshold)
{
  master = merge(em_annotated, de, by.x=0, by.y=0)
  master = master[,-1]
  row.names(master)= master[,"gene_symbols"]
  master$mlog10p = -log10(master[,"p"])
  master$mean_expression = rowMeans(master[,2:10])
  master = master[order(master[,"p"], decreasing=FALSE),]
  master = master[,-(2:10)]
  return(master)
}

#takes a table and returns only the sig genes according to the thresholds provided
de_sig = function(parsed_de_table, p_threshold, fold_threshold)
{
  de_sig_table = subset(parsed_de_table,p.adj < 0.05 & abs(log2fold) > 1.0)
  return(de_sig_table)
}

#plots a volcano plot
plot_volcano = function(de_table, p_threshold, fold_threshold)
{
  de_table_sig = subset(de_table, p.adj <p_threshold & abs(log2fold)>fold_threshold) 
  de_table_sig_up = subset(de_table_sig, p.adj <p_threshold & log2fold>fold_threshold& p.adj!=0)
  de_table_sig_down = subset(de_table_sig, p.adj <p_threshold & log2fold<fold_threshold& p.adj!=0)
  de_table_sig_up = de_table_sig_up[order(de_table_sig_up[,"p"], decreasing=FALSE),]
  de_table_sig_down = de_table_sig_down[order(de_table_sig_down[,"p"], decreasing=FALSE),]
  de_table_sig_up_top10 = de_table_sig_up[1:10,]
  de_table_sig_down_top10 = de_table_sig_down[1:10,]
  de_table_sig_up = subset(de_table_sig, p.adj <p_threshold & log2fold>fold_threshold)
  de_table_sig_down = subset(de_table_sig, p.adj <p_threshold & log2fold<fold_threshold)
  
  ggp = ggplot(de_table, aes(x=log2fold,y=mlog10p)) + 
    geom_point(aes(colour = "black")) +
    geom_point(data = de_table_sig_down, aes(colour="blue")) +
    geom_point(data = de_table_sig_up, aes(colour="red"))+
    labs(title="Volcano Plot", x="log2fold", y="-log10p")+
    my_theme +
    geom_vline(xintercept=-fold_threshold,linetype="dashed", colour = "azure4") +
    geom_vline(xintercept=fold_threshold,linetype="dashed", colour = "azure4") +
    geom_hline(yintercept=-log10(p_threshold),linetype="dashed", colour = "azure4") + 
    xlim(c(-12, 12))+
    ylim(c(0, 300))+
    geom_text_repel(data=de_table_sig_down_top10, aes(label=gene_symbols, colour = "blue"),show.legend = FALSE)+
    geom_text_repel(data=de_table_sig_up_top10, aes(label=gene_symbols, colour = "red"),show.legend = FALSE)+
    scale_colour_manual(values = c("darkslategray", "coral3", "deepskyblue3"),labels=c("non-significant","downregulated","upregulated"))+
    labs(color = "Legend")
  return(ggp)
}

#plots an ma plot
plot_ma = function(de_table, p_threshold, fold_threshold)
{
  de_table_sig = subset(de_table, p.adj <p_threshold & abs(log2fold)>fold_threshold) 
  de_table_sig_up = subset(de_table_sig, p.adj <p_threshold & log2fold>fold_threshold)
  de_table_sig_down = subset(de_table_sig, p.adj <p_threshold & log2fold<fold_threshold)
  de_table_sig_up = de_table_sig_up[order(de_table_sig_up[,"p"], decreasing=FALSE),]
  de_table_sig_down = de_table_sig_down[order(de_table_sig_down[,"p"], decreasing=FALSE),]
  de_table_sig_up_top10 = de_table_sig_up[1:10,]
  de_table_sig_down_top10 = de_table_sig_down[1:10,]
  
  ggp = ggplot(de_table, aes(x=log10(mean_expression),y=log2fold)) + 
    geom_point(aes(colour = "black")) +
    geom_point(data = de_table_sig_down, aes(colour="blue")) +
    geom_point(data = de_table_sig_up, aes(colour="red"))+
    labs(title="MA plot", x="log 10 mean_expression", y="log 2 fold")+
    my_theme +
    geom_hline(yintercept=-2,linetype="dashed", colour = "azure4") + 
    geom_hline(yintercept=2,linetype="dashed", colour = "azure4") + 
    #geom_text_repel(data=de_table_sig_down_top10, aes(label=gene_symbols, colour = "blue"),show.legend = FALSE)+
    #geom_text_repel(data=de_table_sig_up_top10, aes(label=gene_symbols, colour = "red"),show.legend = FALSE)+
    scale_colour_manual(values = c("darkslategray", "coral3", "deepskyblue3"),labels=c("non-significant","downregulated","upregulated"))+
    labs(colour = "Legend")
  return(ggp)
}

#two functions for saving in desired format
save_svg = function(path, ggp)
{
  ggsave(path, plot = ggp)
}
#save_emf = function(path, ggp)
#{
#  emf(path, height=500, width=500)
#  print(ggp)
#  dev.off()
#}

#plots a pca
plot_pca = function(expression_scaled, groups, sample)
{
  EM.nm = as.matrix(sapply(expression_scaled, as.numeric))
  pca = prcomp(t(EM.nm))
  pca_coordinates = data.frame(pca$x)
  
  vars = apply(pca$x, 2, var)
  prop_x = round(vars["PC1"] / sum(vars),4) * 100
  prop_y = round(vars["PC2"] / sum(vars),4) * 100
  x_axis_label = paste("PC1 ", " (",prop_x, "%)",sep="")
  y_axis_label = paste("PC2 ", " (",prop_y, "%)",sep="")
  
  ggp = ggplot(pca_coordinates, aes(x=PC1, y=PC2, colour = groups)) + 
    geom_point(size = 3, shape = 19)+
    scale_color_manual(values=c("aquamarine4", "darkorchid3", "deepskyblue3"))+
    geom_text_repel(size=5, aes(label=sample), show.legend = FALSE)+
    my_theme+
    labs(colour = "Legend")+
    labs(title="PCA", x=x_axis_label, y=y_axis_label)
  return(ggp)
}

#plots the expression density of all samples faceted graph output
plot_exp_density = function(expression_table, columns)
{
  expression_table.m = melt(expression_table)
  expression_table.m$group <- ifelse(grepl("Senes_MtD", expression_table.m$variable), "Senes_MtD", 
                                     ifelse(grepl("Senes", expression_table.m$variable), "Senes", 
                                            ifelse(grepl("Prolif", expression_table.m$variable), "Prolif", "Other")))
  ggp = ggplot(expression_table.m, aes(x = log10(value + 0.001), fill = group)) + 
    geom_density(size = 0.7, alpha = 0.5) +
    my_theme +
    facet_wrap(~variable, ncol = columns) +
    labs(title = "Expression Density Plots", x = "Gene Density", y = "Log 10 Expression") +
    labs(fill = "Legend")+
    scale_fill_manual(values = c("Prolif" = alpha("aquamarine4", 0.3), "Senes" = alpha("darkorchid3", 0.3), "Senes_MtD" = alpha("deepskyblue3", 0.3), "Other" = alpha("gray", 0.3)))
  return(ggp)
}

#makes a single box plot of a gene
make_boxplot = function(gene_name, expression_table, groups)
{
  gene_data = expression_table[gene_name,]
  gene_data = data.frame(t(gene_data))
  gene_data$sample_group = groups
  print(gene_data)
  names(gene_data) = c("expression","sample_group")
  
  ggp = ggplot(gene_data, aes(x=sample_group, y =expression,fill = sample_group)) +
    geom_boxplot(size = 0.5, outlier.size = 0, alpha = 0.5, show.legend=FALSE) +
    geom_jitter(aes(colour =sample_group), width = 0.2, show.legend=FALSE)+
    scale_fill_manual(values = c("Prolif" = "aquamarine4", "Senes" ="darkorchid3", "Senes_MtD" = "deepskyblue3"))+
    my_theme+
    labs(title = "TNF Expression", x = "Sample Groups", y = "Expression")
  return(ggp)
}

#makes a faceted box plot of all genes provided in kist
make_multi_boxplot = function(gene_list, expression_table, groups)
{
  gene_data = expression_table[gene_list,]
  gene_data = data.frame(t(gene_data))
  gene_data$sample_group = groups
  gene_data.m = melt(gene_data, id.vars = "sample_group", variable.name = "Gene", value.name = "Expression")

  ggp = ggplot(gene_data.m, aes(x=sample_group, y =Expression, fill = sample_group)) +
    geom_boxplot(size = 0.5, outlier.size = 0, alpha = 0.5, show.legend=FALSE) +
    scale_fill_manual(values = c("Prolif" = "aquamarine4", "Senes" ="darkorchid3", "Senes_MtD" = "deepskyblue3"))+
    my_theme+
    facet_wrap(~Gene, ncol=10)+
    labs(title = "Multi-Gene Expression", x = "Sample Groups", y = "Expression")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
  return(ggp)
}

#plots a heatmap from scratch, it picks the gene list of the 100 most sig genes in the de table provided and makes the heatmap
plot_heatmap = function(master_sig, em_scaled, cluster_axis, title_given) 
{
  master_sig = master_sig[order(master_sig[,"p"], decreasing=FALSE),]
  sig_genes = row.names(master_sig)
  sig_genes = sig_genes[1:100]
  em_scaled_sig_subset = em_scaled[sig_genes,]
  #part2
  hm.matrix = as.matrix(em_scaled_sig_subset)
  hm.matrix = na.omit(hm.matrix)
  if (cluster_axis == TRUE)
  {
    hm.matrix = t(hm.matrix)
  }
  
  y.dist = Dist(hm.matrix, method="spearman")
  y.cluster = hclust(y.dist, method="average")
  y.dd = as.dendrogram(y.cluster)
  y.dd.reorder = reorder(y.dd,0,FUN="average")
  y.order = order.dendrogram(y.dd.reorder)
  hm.matrix_clustered = hm.matrix[y.order,]
  hm.matrix_clustered = melt(hm.matrix_clustered)
  
  colours = c("deepskyblue3", "cornsilk", "coral3")
  colorRampPalette(colours)(100)
  ggp = ggplot(hm.matrix_clustered, aes(x=Var2, y=Var1, fill =value)) + geom_tile()+ my_theme+
    scale_fill_gradientn(colours = colorRampPalette(colours)(100))+ 
    labs(title= title_given, x = "", y = "Differentially Expressed Genes") +
    theme(axis.text.y = element_blank(), axis.ticks=element_blank(), axis.text.x = element_blank(),legend.title = element_blank(),  panel.border = element_blank())
  return(ggp)
}

#also makes a heatmap but this takes directly the table generated from emscaled and a custom gene list and plots the heatmap
plot_heatmap2 = function(table, title_given)
{
  #clustering and melting the table
  hm.matrix = as.matrix(table)
  hm.matrix = na.omit(hm.matrix)
  y.dist = Dist(hm.matrix, method="spearman")
  y.cluster = hclust(y.dist, method="average")
  y.dd = as.dendrogram(y.cluster)
  y.dd.reorder = reorder(y.dd,0,FUN="average")
  y.order = order.dendrogram(y.dd.reorder)
  hm.matrix_clustered = hm.matrix[y.order,]
  hm.matrix_clustered = melt(hm.matrix_clustered)
  #actually plots the heatmap
  colours = c("deepskyblue3", "cornsilk", "coral3")
  colorRampPalette(colours)(100)
  ggp = ggplot(hm.matrix_clustered, aes(x=Var1, y=Var2, fill =value)) + geom_tile()+ my_theme+
    scale_fill_gradientn(colours = colorRampPalette(colours)(100))+ 
    labs(title= title_given, x = "Differentially Expressed Genes", y = "") +
    theme(axis.text.x = element_blank(), axis.ticks=element_blank(),legend.title = element_blank(),  panel.border = element_blank())
  return(ggp)
}

#plots the rug for a heatmap
plot_rug = function(sample_groups, colour1, colour2)
{
  groups_data = as.matrix(as.numeric(as.factor(sample_groups)))
  groups_data = melt(groups_data)
  colours = c( colour1, colour2)
  ggp = ggplot(groups_data, aes(x = Var1, y = Var2, fill = value)) + 
    geom_tile() + 
    scale_fill_gradientn(colours = colours)+
    geom_tile(linetype="blank")+
    labs( x = "", y = "")+
    theme(legend.position="none", legend.title = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks=element_blank()) 
  return(ggp)
}

#performs pathway analysis and saves all the plots in a directory specified
ora_gsea_plots = function(master,master_sig, p_threshold, fold_threshold, path, em_scaled, groups)
{
  if (!dir.exists(path)) 
  {
    dir.create(path)
  }
  setwd(path)
  colours = c("deepskyblue3", "cornsilk", "coral3")
  colorRampPalette(colours)(100)
  #ora for sig
  sig_genes = row.names(master_sig)
  sig_genes_entrez = bitr(sig_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  ora_results = enrichGO(gene = sig_genes_entrez$ENTREZID, OrgDb = org.Hs.eg.db, readable = T, ont = 
                           "ALL", pvalueCutoff = 0.05, qvalueCutoff = 0.10)
  ggp = barplot(ora_results, showCategory=10)+ 
    scale_fill_gradientn(colours = colorRampPalette(colours)(100))+ 
    theme(axis.text.x = element_text(size = 10))
  save_svg("ora_barplot.svg", ggp)
  
  #ora_for_upregulated
  master_sig_up = subset(master_sig, p.adj <p_threshold & log2fold>fold_threshold)
  sig_genes_up = row.names(master_sig_up)
  sig_genes_entrez_up = bitr(sig_genes_up, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  ora_results = enrichGO(gene = sig_genes_entrez_up$ENTREZID, OrgDb = org.Hs.eg.db, readable = T, ont = 
                           "ALL", pvalueCutoff = 0.05, qvalueCutoff = 0.10)
  ggp = barplot(ora_results, showCategory=10)+ 
    scale_fill_gradientn(colours = colorRampPalette(colours)(100))+ 
    theme(axis.text.x = element_text(size = 8))
  save_svg("ora_barplot_up.svg", ggp)
  
  #ora for downregulated
  master_sig_down = subset(master_sig, p.adj <p_threshold & log2fold<fold_threshold)
  sig_genes_down = row.names(master_sig_down)
  sig_genes_entrez_down = bitr(sig_genes_down, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  ora_results = enrichGO(gene = sig_genes_entrez_down$ENTREZID, OrgDb = org.Hs.eg.db, readable = T, ont = 
                           "ALL", pvalueCutoff = 0.05, qvalueCutoff = 0.10)
  ggp = barplot(ora_results, showCategory=10)+ 
    scale_fill_gradientn(colours = colorRampPalette(colours)(100))+ 
    theme(axis.text.x = element_text(size = 8))
  save_svg("ora_barplot_down.svg", ggp)
  
  #for creating a gsea plot (sorting the values before running gseGO)
  gsea_input = master$log2fold
  names(gsea_input) = row.names(master)
  gsea_input = na.omit(gsea_input)
  gsea_input = sort(gsea_input, decreasing = TRUE)
  gse_results = gseGO(geneList=gsea_input, 
                      ont ="BP", 
                      keyType = "SYMBOL", 
                      nPerm = 10000, 
                      minGSSize = 3, 
                      maxGSSize = 800, 
                      pvalueCutoff = 0.05, 
                      verbose = TRUE, 
                      OrgDb = org.Hs.eg.db, 
                      pAdjustMethod = "none")
  ggp = ridgeplot(gse_results) + 
    scale_fill_gradientn(colours = colorRampPalette(colours)(100))+ 
    theme(axis.text.x = element_text(size = 8))
  save_svg("ridgeplot.svg", ggp)
  
  #ORA enriched genes analysis
  gene_sets = ora_results$geneID
  description = ora_results $Description
  p.adj = ora_results$p.adjust
  ora_results_table = cbind(gene_sets, p.adj)
  ora_results_table = data.frame(ora_results_table)
  rownames(ora_results_table) = description
  write.table(ora_results_table, file="ora_results.tsv", sep="\t")
  enriched_gene_set = as.character(ora_results_table [1,1])
  candidate_genes = unlist(strsplit(enriched_gene_set, "/")) 
  ggp = make_multi_boxplot(candidate_genes[1:10], em_scaled, groups)
  save_svg("ora_boxplot.svg", ggp)
  table_for_heatmap = em_scaled[candidate_genes,]
  data.s = data.frame(t(scale(t(table_for_heatmap))))
  ggp = plot_heatmap2(table_for_heatmap, "ora candidate genes heatmap")
  save_svg("ora_heatmap.svg", ggp)
  
  #GSEA enriched genes analysis
  nes_scores = gse_results$NES
  description = gse_results$Description
  p.adj= gse_results$p.adjust
  gene_sets = gse_results$core_enrichment
  gse_table = cbind(gene_sets, nes_scores, p.adj)
  rownames(gse_table) = description
  write.table(gse_table, file="gse_results.tsv", sep="\t")
  enriched_gene_set = as.character(gse_table[1,1])
  candidate_genes = unlist(strsplit(enriched_gene_set, "/")) 
  ggp = make_multi_boxplot(candidate_genes[1:10], em_scaled, groups)
  save_svg("gsea_boxplot.svg", ggp)
  table_for_heatmap = em_scaled[candidate_genes,]
  data.s = data.frame(t(scale(t(table_for_heatmap))))
  ggp = plot_heatmap2(table_for_heatmap, "gsea candidate genes heatmap")
  save_svg("gsea_heatmap.svg", ggp)
  
  #STRING Network (using gsea most enriched only)
  candidate_genes_table = data.frame(candidate_genes)
  names(candidate_genes_table) = "gene"
  string_db = STRINGdb$new( version="11.5", species=9606, score_threshold=200, 
                            network_type="full", input_directory="")
  string_mapped = string_db$map(candidate_genes_table, "gene", removeUnmappedRows = TRUE )
  ggp = string_db$plot_network(string_mapped)
  save_svg("gsea_string.svg", ggp)
}

#makes a metagene multi box plot using the signature identified from the heatmap of all sig genes
make_metagene_boxplot = function(signature, em_annotated, groups) 
{
  em_symbols = em_annotated[,2:10]
  
  #scaling again with gene_id in rownames(the main em_scaled in the code has symbols in rownames)
  em_scaled = t(em_symbols)
  em_scaled = data.frame(em_scaled)
  em_scaled = scale(em_scaled)
  em_scaled = t(em_scaled)
  em_scaled = data.frame(em_scaled)
  em_scaled = na.omit(em_scaled)
  
  #metagene
  intermediate_table = em_scaled[signature,]
  intermediate_table = na.omit(intermediate_table)
  metagene = data.frame(colMeans(intermediate_table))
  names(metagene) = "value"
  metagene$group = groups
  metagene$group = factor(metagene$group, levels=c("Prolif", "Senes", "Senes_MtD"))
  # plot
  ggp = ggplot(metagene, aes(x=group, y= value, colour = group, fill= group)) + 
    geom_boxplot()+
    labs(fill = "Legend", colour = "Legend")+
    my_theme+
    labs(title = "Metagene Plot", x= "Groups", y = "Mean Z-Score")
  return(ggp)
}

