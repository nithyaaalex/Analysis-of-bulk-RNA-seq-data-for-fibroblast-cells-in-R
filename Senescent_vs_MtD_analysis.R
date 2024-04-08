#### LOADING LIBRARIES ####
library(ggplot2)
library(ggrepel)
library(reshape2)
library(amap)
library(devEMF)
library(svglite)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggridges)
library(STRINGdb)
source("C:/Users/2933044s/Desktop/assessment/functions_for_analysis.R")



#### LOAD DATA ####
em = read.table("C:/Users/2933044s/Desktop/assessment/EM.csv", header=TRUE, row.names = 1, sep= "\t")
ss = read.table("C:/Users/2933044s/Desktop/assessment/sample_sheet.csv", header=TRUE, sep= "\t")
annotations = read.table("C:/Users/2933044s/Desktop/assessment/Human_Background_GRCh38.p13.csv", header=TRUE, row.names = 1, sep= "\t")
names(annotations) =  c("gene_symbols","chromosome","start","stop","bio_type")
de_mtd_vs_prolif = read.table("C:/Users/2933044s/Desktop/assessment/DE_Senes_MtD_vs_Prolif.csv", header=TRUE, row.names = 1, sep= "\t")
de_mtd_vs_senes = read.table("C:/Users/2933044s/Desktop/assessment/DE_Senes_MtD_vs_Senes.csv", header=TRUE, row.names = 1, sep= "\t")
de_senes_vs_prolif = read.table("C:/Users/2933044s/Desktop/assessment/DE_Senes_vs_Prolif.csv", header=TRUE, row.names = 1, sep= "\t")
##declaring global thresholds##
p_threshold = 0.001
fold_threshold = 2


#### PARSE DATA ####

##creating master table for complex analysis first##
#adding sig columns to all de tables -  doing this before the parse_de function because need it for master too
de_mtd_vs_prolif$sig = factor(de_mtd_vs_prolif$p.adj < p_threshold & abs(de_mtd_vs_prolif$log2fold) > fold_threshold)
de_mtd_vs_senes$sig = factor(de_mtd_vs_senes$p.adj < p_threshold & abs(de_mtd_vs_senes$log2fold) > fold_threshold)
de_senes_vs_prolif$sig = factor(de_senes_vs_prolif$p.adj < p_threshold & abs(de_senes_vs_prolif$log2fold) > fold_threshold)
#merging them
master_temp = merge(de_mtd_vs_prolif, de_mtd_vs_senes, by.x=0, by.y=0, suffixes=c(".mvp",".mvs"))
row.names(master_temp) = master_temp$"Row.names"
master_temp = master_temp [,-1]
master = merge(master_temp, de_senes_vs_prolif,by.x=0, by.y=0, suffixes = c("",".svp"))
#the no suffix values are .svp - the suffix isn't getting added#
row.names(master) = master$"Row.names" 
master = master[, -1]


#creating em_annotated to merge with the de tables#
em_annotated = merge(em, annotations, by.x=0, by.y=0)
names(em_annotated)[1] = "gene_id"
row.names(em_annotated) = em_annotated[,"gene_id"]

#using parse_de function to parse all 3 de files and annotate them#
de_mtd_vs_prolif = parse_de(de_mtd_vs_prolif, em_annotated, p_threshold, fold_threshold)
de_mtd_vs_senes = parse_de(de_mtd_vs_senes, em_annotated, p_threshold, fold_threshold)
de_senes_vs_prolif = parse_de(de_senes_vs_prolif, em_annotated, p_threshold, fold_threshold)

#using de_sig to only save the significant genes#
de_mtd_vs_prolif_sig = de_sig(de_mtd_vs_prolif, p_threshold, fold_threshold)
de_mtd_vs_senes_sig = de_sig(de_mtd_vs_senes, p_threshold, fold_threshold)
de_senes_vs_prolif_sig = de_sig(de_senes_vs_prolif, p_threshold, fold_threshold)

#creating expression symbols table and naming with gene_symbols#
em_symbols = em_annotated[,2:11]
row.names(em_symbols) = em_symbols[,"gene_symbols"]
em_symbols = em_symbols[,-10]
row.names(em_annotated) = em_annotated[,"gene_symbols"]

#creating scaled table for further use#
em_scaled = t(em_symbols)
em_scaled = data.frame(em_scaled)
em_scaled = scale(em_scaled)
em_scaled = t(em_scaled)
em_scaled = data.frame(em_scaled)
em_scaled = na.omit(em_scaled)

#subsetting ss for rugs#
ss_prolif_vs_senes = subset(ss[1:6,])
ss_mtd_vs_prolif = subset(ss[c(1:3,7:9),])
ss_mtd_vs_senes = subset(ss[4:9,])



#### THEME ####
my_theme = theme(
  text = element_text(family = "sans"),
  plot.title = element_text(size=23),
  axis.text.x = element_text(size=13),
  axis.text.y = element_text(size=13),
  axis.title.x = element_text(size=19),
  axis.title.y = element_text(size=19),
  legend.text = element_text(size=13),
  legend.title = element_text(size = 17, face = "bold"),
  legend.position = "right",
  legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
  panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
  panel.background = element_rect(fill = "white"),
  panel.grid = element_blank(),
  panel.spacing = unit(1, "lines"),
  strip.text.x = element_text(size = 16, family="Arial", face="bold", vjust=1)
)



#### PLOTS ####
#my colours for samples "aquamarine4", "darkorchid3", "deepskyblue3"#
#my colours for non-sig, up, down "darkslategray", "coral", "cornflowerblue"#

###generic plots###
ggp = plot_pca(em_scaled, ss$SAMPLE_GROUP, ss$SAMPLE)
save_svg("C:/Users/2933044s/Desktop/assessment/pca.svg", ggp)
ggp = plot_exp_density(em, 3)
save_svg("C:/Users/2933044s/Desktop/assessment/exp_density.svg", ggp)


###ma plots###
#ggp = plot_ma(de_mtd_vs_prolif, p_threshold, fold_threshold)
#ggp = plot_ma(de_mtd_vs_senes , p_threshold, fold_threshold)
ggp = plot_ma(de_senes_vs_prolif, p_threshold, fold_threshold)
save_svg("C:/Users/2933044s/Desktop/assessment/ma(senes_prolif).svg", ggp)


###volcano plots###
ggp = plot_volcano(de_mtd_vs_prolif, p_threshold, fold_threshold)
ggp = plot_volcano(de_mtd_vs_senes , p_threshold, fold_threshold)
ggp = plot_volcano(de_senes_vs_prolif, p_threshold, fold_threshold)
save_svg("C:/Users/2933044s/Desktop/assessment/volcano(senes_prolif).svg", ggp)


##plot for comparing total number of differentially expressed genes##
values = c(nrow(subset(de_mtd_vs_prolif_sig, sig == TRUE & log2fold>fold_threshold)),
nrow(subset(de_mtd_vs_prolif_sig, sig == TRUE& log2fold<fold_threshold)),
nrow(subset(de_mtd_vs_senes_sig, sig == TRUE& log2fold>fold_threshold)),
nrow(subset(de_mtd_vs_senes_sig, sig  == TRUE & log2fold<fold_threshold)),
nrow(subset(de_senes_vs_prolif_sig, sig == TRUE& log2fold>fold_threshold)),
nrow(subset(de_senes_vs_prolif_sig, sig == TRUE& log2fold<fold_threshold)))
labels = c("prolif_up", "prolif_down", "mtd_up", "mtd_down", "senes_up", "senes_down" )
groups = c("up", "down", "up", "down", "up", "down" )
temp_table = data.frame(Column1 = values, Column2 = labels, Column3 = groups)
ggp = ggplot(temp_table, aes(x = labels, y = values, fill = groups)) + 
  geom_bar(stat = "identity", position = "dodge") +
  xlab("Levels") +
  ylab("Count") +
  labs(fill = "Legend") +
  my_theme+
  ggtitle("Bar Plot of Upregulated and Downregulated genes")+
  scale_fill_manual(values = c("deepskyblue3", "coral3"))+
  theme(plot.title = element_text(size=17))
save_svg("C:/Users/2933044s/Desktop/assessment/overall_number_bar_plot.svg", ggp)
ggp

###multi-box-plots###
#gene_list = row.names(de_mtd_vs_prolif[1:10,])
#gene_list = row.names(de_mtd_vs_senes [1:10,])
#all were the same genes selected - save using export
gene_list = row.names(de_senes_vs_prolif[1:10,])
ggp = make_multi_boxplot(gene_list, em_scaled, ss$SAMPLE_GROUP)
ggp


###heatmaps###
ggp = plot_heatmap(de_mtd_vs_prolif, em_scaled[,c(1:3,7:9)], FALSE, "Prolif vs Senes_MtD")
save_svg("C:/Users/2933044s/Desktop/assessment/heatmap(mtd_prolif).svg", ggp)
ggp = plot_heatmap(de_mtd_vs_senes, em_scaled[,4:9], FALSE, "Senes vs Senes_MtD")
save_svg("C:/Users/2933044s/Desktop/assessment/heatmap(mtd_senes).svg", ggp)
ggp = plot_heatmap(de_senes_vs_prolif, em_scaled[1:6], FALSE, "Prolif vs Senes")
save_svg("C:/Users/2933044s/Desktop/assessment/heatmap(senes_prolif).svg", ggp)


##rugs##
mtd = "deepskyblue3" #setting colours for functions
prolif = "aquamarine4"
senes = "darkorchid3"
ggp = plot_rug(ss_mtd_vs_prolif$SAMPLE_GROUP, prolif, mtd)
save_svg("C:/Users/2933044s/Desktop/assessment/rug(mtd_prolif).svg", ggp)
ggp = plot_rug(ss_mtd_vs_senes$SAMPLE_GROUP, senes, mtd)
save_svg("C:/Users/2933044s/Desktop/assessment/rug(mtd_senes).svg", ggp)
ggp = plot_rug(ss_prolif_vs_senes$SAMPLE_GROUP, prolif, senes)
save_svg("C:/Users/2933044s/Desktop/assessment/rug(senes_prolif).svg", ggp)


###pathway plots###
ora_gsea_plots(de_mtd_vs_prolif, de_mtd_vs_prolif_sig, p_threshold, fold_threshold, "C:/Users/2933044s/Desktop/assessment/mtd_vs_prolif",em_scaled, ss$SAMPLE_GROUP)
ora_gsea_plots(de_mtd_vs_senes, de_mtd_vs_senes_sig, p_threshold, fold_threshold, "C:/Users/2933044s/Desktop/assessment/mtd_vs_senes",em_scaled, ss$SAMPLE_GROUP)
ora_gsea_plots(de_senes_vs_prolif, de_senes_vs_prolif_sig, p_threshold, fold_threshold, "C:/Users/2933044s/Desktop/assessment/prolif_vs_senes",em_scaled, ss$SAMPLE_GROUP)


###complex plots###
##foldvsfold##
cor = cor.test(master$log2fold.mvs, master$log2fold.mvp)
ggp = ggplot(master, aes(x = log2fold.mvs, y = log2fold.mvp)) +
      geom_point() +
      labs(title = paste("SenesVsSenesMtd and SenesMtDvsProlif",cor), x = "Log2Fold", y = "Log2Fold")
save_svg("C:/Users/2933044s/Desktop/assessment/SenesVsSenesMtd and SenesMtDvsProlif.svg", ggp)
cor = cor.test(master$log2fold, master$log2fold.mvp)
ggp = ggplot(master, aes(x = log2fold, y = log2fold.mvp)) +
  geom_point() +
  labs(title = paste("SenesVsProlif and SenesMtDvsProlif",cor), x = "Log2Fold", y = "Log2Fold")
save_svg("C:/Users/2933044s/Desktop/assessment/SenesVsProlif and SenesMtDvsProlif.svg", ggp)
cor = cor.test(master$log2fold.mvs, master$log2fold)
ggp = ggplot(master, aes(x = log2fold.mvs, y = log2fold)) +
  geom_point() +
  labs(title = paste("SenesVsSenesMtd and SenesvsProlif",cor), x = "Log2Fold", y = "Log2Fold")
save_svg("C:/Users/2933044s/Desktop/assessment/SenesVsSenesMtd and SenesvsProlif.svg", ggp)


##heatmap (of ALL sig genes)##

#gets all the significant genes and makes table for it
all_sig = subset(master, sig.mvp == TRUE | sig.mvs == TRUE | sig == TRUE )
mega_gene_list = row.names(all_sig)
sig_all_expression = em_scaled[mega_gene_list,]
data.s = data.frame(t(scale(t(sig_all_expression))))
ggp = plot_heatmap2(data.s, "Heatmap of ALL significant genes")
save_svg("C:/Users/2933044s/Desktop/assessment/ALLheatmap.svg", ggp)

##differential signatures##
#signature_1 = rownames(subset(master, (sig.mvp == FALSE) & (sig.mvs == TRUE & log2fold.mvs < 0)& (sig == TRUE & log2fold >0)))
#signature_2 = rownames(subset(master, (sig.mvp == TRUE & log2fold.mvp < 0) & (sig.mvs == FALSE)& (sig == TRUE & log2fold <0)))
#signature_3 = rownames(subset(master, (sig.mvp == TRUE & log2fold.mvp > 0) & (sig.mvs == TRUE & log2fold.mvs > 0)& (sig == FALSE)))
                     
                                         
##metagene plots##
#ggp = make_metagene_boxplot(signature_1, em_annotated,ss$SAMPLE_GROUP)
#ggp = make_metagene_boxplot(signature_2, em_annotated,ss$SAMPLE_GROUP)
#ggp = make_metagene_boxplot(signature_3, em_annotated,ss$SAMPLE_GROUP)
