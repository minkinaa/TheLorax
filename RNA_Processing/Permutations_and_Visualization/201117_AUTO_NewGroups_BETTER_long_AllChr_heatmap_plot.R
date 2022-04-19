##This is a version of 191113_Tree29_make_group_by_gene_files.R, updated with new (final!!) lineage groups.

###This one makes the long plot but starting with saved table

### Input file "All_NonLog_Normalized_Scaled_gene_counts" is the raw counts file (available in GEO) scaled so that every cell has 10,000 counts

library(dplyr)
#library(Seurat)
library(ggplot2)
library(ggrepel)
library(reshape2)
library(robustbase) #colMedians
library("RColorBrewer")
library(gridExtra)
library(matrixStats)

#Output_dir = "/Users/anna/Documents/UW_GS/Shendure_Lab/Projects/Many_targets_ideas/CROPt/190722_Tree29_other7plates/TARG_plates2thru8/analysis/all_plates/191106_New_Groups_from_Run13_take2/Reordered_group_plots/"
#Output_dir = "/Users/anna/Documents/UW_GS/Shendure_Lab/Projects/Many_targets_ideas/CROPt/190722_Tree29_other7plates/TARG_plates2thru8/200623_NewPipelineAnalysis/201117_All_RNA_Analysis/201210_LongPlots_Rev_Order/"

args = commandArgs(trailingOnly=TRUE)
Output_dir = args[1]

gene_metadata_with_names = read.table("/Users/anna/Documents/UW_GS/Shendure_Lab/Projects/Many_targets_ideas/CROPt/CROPt_pipeline_files/Commonly_sourced_files/unique_gencode.v29.gene_ids_gene_names_PLUS.bed", header = TRUE)
gene_metadata_with_names$V2 = gene_metadata_with_names$undated_start_pos
chr_pos = read.table("/Users/anna/Documents/UW_GS/Shendure_Lab/Projects/Many_targets_ideas/CROPt/CROPt_pipeline_files/Commonly_sourced_files/chrom_length_midpoint_etc", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

#group_file = as.data.frame(read.table("/Users/anna/Documents/UW_GS/Shendure_Lab/Projects/Many_targets_ideas/CROPt/190722_Tree29_other7plates/TARG_plates2thru8/200623_NewPipelineAnalysis/201116_PostTreeBuild_AddAll_Cells/All_Cells_W_LinGroups_And_Rep_Cells_AND_NEW_Rep_Cells", header = TRUE, sep = "\t"), drop = FALSE)
group_file_name = args[2]
group_file =  as.data.frame(read.table(group_file_name, sep = "\t", stringsAsFactors = FALSE))
reorder = TRUE
high_genes = "keep" # "remove" is other option

#new_order = c(8,7,4,1,9,3,2,10,5,6,30,31,29,28,24,26,25,23,41,42,56,20,21,22,18,57,17,15,16,19,14,11,13,12,44,45,40,50,55,52,53,54,51,34,36,35,27,46,39,38,43,48,37,49,32,33,47)
new_order = c(1:max(group_file$V2))
new_order_char = paste0("grp", new_order)

#I *think* this is the file to be used in the next step 
#test_file = read.table("/Users/anna/Documents/UW_GS/Shendure_Lab/Projects/Many_targets_ideas/CROPt/190722_Tree29_other7plates/TARG_plates2thru8/analysis/all_plates/191106_New_Groups_from_Run13_take2/All_NonLog_Normalized_Scaled_gene_counts", sep = "\t", row.names = 1)
#test_file = test_file[match(rownames(group_file), rownames(test_file)),]

print("Reading Large Table")

#NonLog_Nor_counts = read.table("/Users/anna/Documents/UW_GS/Shendure_Lab/Projects/Many_targets_ideas/CROPt/190722_Tree29_other7plates/TARG_plates2thru8/analysis/all_plates/191106_New_Groups_from_Run13_take2/All_NonLog_Normalized_Scaled_gene_counts", sep = "\t", row.names = 1)
NonLog_Nor_counts = read.table("/Users/anna/Documents/UW_GS/Shendure_Lab/Projects/Many_targets_ideas/CROPt/190722_Tree29_other7plates/TARG_plates2thru8/200623_NewPipelineAnalysis/201117_All_RNA_Analysis/201117_New_Cell_by_Gene_Scaled_etc_Files/All_NonLog_Normalized_Scaled_gene_counts", sep = "\t", row.names = 1)
NonLog_Nor_counts = NonLog_Nor_counts[match(rownames(group_file), rownames(NonLog_Nor_counts)), which(colnames(NonLog_Nor_counts) %in% as.character(gene_metadata_with_names$V7))]
###The above code gets rid of any cells which got tossed!

#test = NonLog_Nor_counts[which(rownames(NonLog_Nor_counts) %in% rownames(group_file)),]

#all_orig_cells = read.table("/Users/anna/Documents/UW_GS/Shendure_Lab/Projects/Many_targets_ideas/CROPt/190722_Tree29_other7plates/TARG_plates2thru8/analysis/all_plates/all_plates_cellranger_files_2/barcodes.tsv", sep = "\t", stringsAsFactors = FALSE)
#all_orig_cells$V2 = all_orig_cells$V1
#all_orig_cells = all_orig_cells[which(all_orig_cells$V1 %in% rownames(group_file)),]

NonLog_Nor_counts$Temp_col = group_file$V2

print("Calculating Means, Etc.")

all_genes_on_chr = gene_metadata_with_names[,c(1,2,3,7)]
all_genes_on_chr_in_data = all_genes_on_chr[which(all_genes_on_chr$V7 %in% colnames(NonLog_Nor_counts)),]
all_genes_on_chr_in_data = all_genes_on_chr_in_data[match(as.character(all_genes_on_chr_in_data$V7), colnames(NonLog_Nor_counts[,1:ncol(NonLog_Nor_counts)-1])),]
all_genes_on_chr_in_data$mean_exp_in_data = colMeans(NonLog_Nor_counts[,1:ncol(NonLog_Nor_counts)-1])
all_genes_on_chr_in_data$median_exp_in_data = colMedians(as.matrix(NonLog_Nor_counts[,1:ncol(NonLog_Nor_counts)-1]))
all_genes_on_chr_in_data$var_exp_in_data = colVars(as.matrix(NonLog_Nor_counts[,1:ncol(NonLog_Nor_counts)-1]))

plot(all_genes_on_chr_in_data[which(all_genes_on_chr_in_data$mean_exp_in_data < 25),]$mean_exp_in_data, all_genes_on_chr_in_data[which(all_genes_on_chr_in_data$mean_exp_in_data < 25),]$median_exp_in_data)
plot(all_genes_on_chr_in_data[which(all_genes_on_chr_in_data$mean_exp_in_data < 15),]$mean_exp_in_data, all_genes_on_chr_in_data[which(all_genes_on_chr_in_data$mean_exp_in_data < 15),]$var_exp_in_data)
plot(all_genes_on_chr_in_data[which(all_genes_on_chr_in_data$mean_exp_in_data < 15),]$mean_exp_in_data, (all_genes_on_chr_in_data[which(all_genes_on_chr_in_data$mean_exp_in_data < 15),]$var_exp_in_data/all_genes_on_chr_in_data[which(all_genes_on_chr_in_data$mean_exp_in_data < 15),]$mean_exp_in_data))

all_group_means = data.frame()
all_genes_plus_group_means = all_genes_on_chr_in_data[,1:4]

for (i in 1:max(NonLog_Nor_counts$Temp_col))
{
  group_num = i
  print (group_num)
  NonLog_Nor_counts_subgroup = NonLog_Nor_counts[which(NonLog_Nor_counts$Temp_col == i),]
  all_genes_plus_group_means$temp_means = colMeans(NonLog_Nor_counts_subgroup[,1:ncol(NonLog_Nor_counts_subgroup)-1])
  if (i > 1) {
    all_group_means = cbind(all_group_means, temp_name = as.data.frame(all_genes_plus_group_means[,5]))
  } else
  {
    all_group_means = as.data.frame(all_genes_plus_group_means[,5], drop = FALSE)
    colnames(all_group_means) = paste0("grp",i)
  }
}
colnames(all_group_means) = paste0("grp", c(1:ncol(all_group_means)))
rownames(all_group_means) = all_genes_on_chr_in_data$V7

write.table(all_group_means, paste0(Output_dir,"NonLog_All_groups_by_gene_means_Updated"), quote = FALSE, sep = "\t")

if (reorder == TRUE)
{
  all_group_means = all_group_means[,match(new_order_char,colnames(all_group_means))]
}

all_group_means$max_mean = apply(all_group_means, 1, FUN = max)
#all_group_means$median = apply(all_group_means[,1:57], 1, FUN = median)

mean_cutoff = .5

all_group_means_cutoff = all_group_means[which(all_group_means$max_mean > mean_cutoff),]
all_genes_on_chr_in_data_cutoff = all_genes_on_chr_in_data[which(all_genes_on_chr_in_data$V7 %in% rownames(all_group_means_cutoff)),]

all_group_means_cutoff_div_mean = all_group_means_cutoff/all_genes_on_chr_in_data_cutoff$mean_exp_in_data
all_group_means_cutoff_div_mean_log = log(all_group_means_cutoff_div_mean, base = 2)
all_group_means_cutoff_div_mean_log = all_group_means_cutoff_div_mean_log[,1:ncol(all_group_means_cutoff_div_mean_log)-1]

all_group_means_cutoff_div_median = all_group_means_cutoff/all_genes_on_chr_in_data_cutoff$median_exp_in_data
all_group_means_cutoff_div_median_log = log(all_group_means_cutoff_div_median, base = 2)
all_group_means_cutoff_div_median_log = all_group_means_cutoff_div_median_log[,1:ncol(all_group_means_cutoff_div_median_log)-1]

all_group_means_cutoff_MINUS_mean = all_group_means_cutoff - all_genes_on_chr_in_data_cutoff$mean_exp_in_data
all_group_means_cutoff_MINUS_mean = all_group_means_cutoff_MINUS_mean[,1:ncol(all_group_means_cutoff_MINUS_mean)-1]

all_group_means_cutoff_MINUS_mean_scaled = t(scale(t(all_group_means_cutoff_MINUS_mean[1:(ncol(all_group_means_cutoff_MINUS_mean))])))
all_group_means_cutoff_MINUS_mean_scaled = all_group_means_cutoff_MINUS_mean_scaled[,1:ncol(all_group_means_cutoff_MINUS_mean_scaled)]

mean_or_median = "mean"
temp_min = -1.5
temp_max = 1.5

if (mean_or_median == "mean")
{
  temp_matrix = as.matrix(all_group_means_cutoff_div_mean_log)
} else if (mean_or_median == "median") {
  temp_matrix = as.matrix(all_group_means_cutoff_div_median_log)
} else if (mean_or_median == "minus_mean") {
  temp_matrix = as.matrix(all_group_means_cutoff_MINUS_mean)
} else if (mean_or_median == "minus_mean_scaled") {
  temp_matrix = as.matrix(all_group_means_cutoff_MINUS_mean_scaled)
}

#write.table(all_group_means_cutoff_MINUS_mean_scaled, paste0(Output_dir, "Table_MinusMinScaled_NonLog_-4to4_cutoff1"), sep = "\t", quote = FALSE)

temp_matrix_melt = melt(temp_matrix)
temp_matrix_melt$width = 1
temp_matrix_melt$height = 1
hist(temp_matrix_melt$value, breaks = 100)
hist(temp_matrix_melt[which(temp_matrix_melt$Var1 == "TMEM135"),]$value, breaks = 100)

bad_genes_temp = temp_matrix_melt[which(temp_matrix_melt$value < temp_min | temp_matrix_melt$value > temp_max),]

###### playing with new stuff #######

write.table(temp_matrix, paste0(Output_dir, "Matrix_File_W_Means"), sep = "\t", quote = FALSE)

# tree_coordinates_file = read.table("/Users/anna/Documents/UW_GS/Shendure_Lab/Projects/Many_targets_ideas/CROPt/190722_Tree29_other7plates/TARG_plates2thru8/200623_NewPipelineAnalysis/201116_PostTreeBuild_AddAll_Cells/Consensus_Plot_New_Counts/Tree_Coordinates_File", sep = "\t", stringsAsFactors = FALSE, header = TRUE)
# Output_dir = "/Users/anna/Documents/UW_GS/Shendure_Lab/Projects/Many_targets_ideas/CROPt/190722_Tree29_other7plates/TARG_plates2thru8/200623_NewPipelineAnalysis/201117_All_RNA_Analysis/201210_PlayingW_Tree_Plus_Long_Plot/"
# 
# scale_x = .25
# tree_coordinates_file$x1 = scale_x*tree_coordinates_file$x1
# tree_coordinates_file$x2 = scale_x*tree_coordinates_file$x2
# 
# gene_nums = as.data.frame(cbind(rownames(temp_matrix), (1:nrow(temp_matrix))+max(tree_coordinates_file$x2)))
# colnames(gene_nums) = c("Gene", "Pos")
# 
# groups_and_nums = as.data.frame(cbind(colnames(temp_matrix), (ncol(temp_matrix):1)-.5))
# colnames(groups_and_nums) = c("Group", "Y_Pos")
# 
# new_temp_matrix = temp_matrix
# rownames(new_temp_matrix) = gene_nums$Pos
# colnames(new_temp_matrix) = groups_and_nums$Y_Pos
# new_temp_matrix_melt = melt(new_temp_matrix)
# 
# new_temp_matrix_melt$width = 1
# new_temp_matrix_melt$height = 1
# 
# new_gene_subset = gene_nums[which(gene_nums$Gene %in% gene_subset$V7),]
# new_temp_matrix_melt_subset = new_temp_matrix_melt[which(as.character(new_temp_matrix_melt$Var1) %in% new_gene_subset$Pos),]
# 
# min_heatmap_x = min(new_temp_matrix_melt_subset$Var1)
# max_heatmap_x = max(tree_coordinates_file$x2)
# 
# diff = min_heatmap_x-max_heatmap_x
# 
# new_temp_matrix_melt_subset$Var1 = new_temp_matrix_melt_subset$Var1 - diff
# 
# 
# TEST_plot_tbl_ordered = ggplot(new_temp_matrix_melt_subset, aes(x = Var1, y = Var2)) +
#   geom_tile(aes(fill=value), width = new_temp_matrix_melt_subset$width, height = new_temp_matrix_melt_subset$height) + 
#   scale_fill_gradientn(colours=colors_temp, breaks = c(temp_min,-.05,.05,temp_max), limits=c(temp_min,temp_max)) +
#   theme_minimal() + geom_segment(data = tree_coordinates_file, aes(x = x1, y = y1, xend = x2, yend = y2, colour = "segment")) +
#   theme(legend.position="none") +
#   theme(axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x=element_blank()) +
#   labs(title = plot_title)
# 
# pdf(paste0(Output_dir,paste0("Long_plot_NonLog_",mean_or_median,"_",temp_min,"to",temp_max,"_cutoff_",mean_cutoff,"_",high_genes,"_", plot_title,".pdf"), sep = ""), width = (num_genes/40)+.5, height = 10)
# print(TEST_plot_tbl_ordered )
# dev.off()


###### playing with new stuff ABOVE #######


if (high_genes == "remove")
{
  temp_matrix_melt = temp_matrix_melt[which(!(as.character(temp_matrix_melt$Var1) %in% as.character(bad_genes_temp$Var1))),]
} else if (high_genes == "keep") {
  if(nrow(temp_matrix_melt[which(temp_matrix_melt$value > temp_max),]) > 0){
    temp_matrix_melt[which(temp_matrix_melt$value > temp_max),]$value = temp_max
  }
  if(nrow(temp_matrix_melt[which(temp_matrix_melt$value < temp_min),]) > 0){
    temp_matrix_melt[which(temp_matrix_melt$value < temp_min),]$value = temp_min
  }
}

#temp_matrix_melt = temp_matrix_melt[which(!(as.character(temp_matrix_melt$Var1) %in% as.character(bad_genes_temp$Var1))),]
length(unique(temp_matrix_melt$Var1))
#temp_matrix_melt = temp_matrix_melt[which((as.character(temp_matrix_melt$Var1) %in% as.character(bad_genes_temp$Var1))),]

colors_temp = rev(brewer.pal(11,"RdBu"))

List_of_plots = list()
plot_lengths_vect = vector()
chrom_names = unique(as.character(all_genes_on_chr_in_data$V1))

for (i in 1:length(chrom_names))
{
  print (i)
  chrom = chrom_names[i]
  plot_title = chrom
  gene_subset = all_genes_on_chr_in_data[which(all_genes_on_chr_in_data$V1 == chrom),]
  temp_matrix_melt_subset = temp_matrix_melt[which(as.character(temp_matrix_melt$Var1) %in% gene_subset$V7),]
  plot_lengths_vect = c(plot_lengths_vect, rep(i,length(unique(temp_matrix_melt_subset$Var1))))
  if (i < 23)
  {
    plot_tbl_ordered = ggplot(temp_matrix_melt_subset, aes(x = Var1, y = Var2)) +
      geom_tile(aes(fill=value), width = temp_matrix_melt_subset$width, height = temp_matrix_melt_subset$height) + 
      scale_fill_gradientn(colours=colors_temp, breaks = c(temp_min,-.05,.05,temp_max), limits=c(temp_min,temp_max)) +
      theme_minimal() + ylim(rev(levels(temp_matrix_melt_subset$Var2))) +
      theme(legend.position="none") +
      theme(axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x=element_blank()) +
      labs(title = plot_title)
    #print(plot_tbl_ordered)
    List_of_plots[[i]] = plot_tbl_ordered
  } else #keep legend
  {
    plot_tbl_ordered = ggplot(temp_matrix_melt_subset, aes(x = Var1, y = Var2)) +
      geom_tile(aes(fill=value), width = temp_matrix_melt_subset$width, height = temp_matrix_melt_subset$height) + 
      scale_fill_gradientn(colours=colors_temp, breaks = c(temp_min,-.05,.05,temp_max), limits=c(temp_min,temp_max)) +
      theme_minimal() + ylim(rev(levels(temp_matrix_melt_subset$Var2))) +
      theme(axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x=element_blank()) +
      labs(title = plot_title)
    #print(plot_tbl_ordered)
    List_of_plots[[i]] = plot_tbl_ordered
  }
  #generate and save individual chromosome plot
  plot_individual = ggplot(temp_matrix_melt_subset, aes(x = Var1, y = Var2)) +
    geom_tile(aes(fill=value), width = temp_matrix_melt_subset$width, height = temp_matrix_melt_subset$height) + 
    scale_fill_gradientn(colours=colors_temp, breaks = c(temp_min,-.05,.05,temp_max), limits=c(temp_min,temp_max)) +
    theme_minimal() +
    theme(legend.position="none") +
    theme(axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(title = plot_title)
  num_genes = length(unique(temp_matrix_melt_subset$Var1))
  pdf(paste0(Output_dir,paste0("Long_plot_NonLog_",mean_or_median,"_",temp_min,"to",temp_max,"_cutoff_",mean_cutoff,"_",high_genes,"_", plot_title,".pdf"), sep = ""), width = (num_genes/5)+.5, height = 10)
  print(plot_individual)
  dev.off()
}
lay = t(matrix(plot_lengths_vect))
pdf(paste0(Output_dir,paste0("Long_plot_NonLog_",mean_or_median,"_",temp_min,"to",temp_max,"_cutoff_",mean_cutoff,"_", high_genes,".pdf"), sep = ""), width = (150*(dim(lay)[2])/5000), height = (10))
#pdf(paste0(Output_dir,paste0("Opp_Long_plot_NonLog_",mean_or_median,"_",temp_min,"to",temp_max,"_cutoff_",mean_cutoff,".pdf"), sep = ""), width = 40, height = 10)
print(grid.arrange(grobs = List_of_plots, layout_matrix = lay, ncol = length(chrom_names)))
dev.off()




        

