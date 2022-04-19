source("/Users/anna/Documents/UW_GS/Shendure_Lab/Projects/Many_targets_ideas/CROPt/CROPt_pipeline_files/Scripts/R_FUNCTIONS_general_header_file.R")
library(reshape2)
library(RColorBrewer)
library(ggplot2)
library(gridExtra)
options(stringsAsFactors = FALSE)

Output_dir = "/Users/anna/Documents/UW_GS/Shendure_Lab/Projects/Many_targets_ideas/CROPt/191126_targATAC_full_1/210405_Tree42_ATAC_Analysis/bin5000upstreamPlusGeneBody/LongPlots/"
bin_size = "bin5000upstreamPlusGeneBody"
#bin_size = "1000000"
scaled = "NOT_SCALED"
prefix = paste0(bin_size, "bins_", scaled)

binned_by_gene = FALSE
just_gene_name_no_bins = TRUE
only_keep_set_list_of_genes = FALSE

encode_bins = FALSE

#ATAC_bins_tbl_name = "/Users/anna/Documents/UW_GS/Shendure_Lab/Projects/Many_targets_ideas/CROPt/191126_targATAC_full_1/210405_Tree42_ATAC_Analysis/bin2000upst0downst/AllPlates_upstr2000_dwnstr0_JustReadsInPeaks.txt"  
ATAC_bins_tbl_name = "/Users/anna/Documents/UW_GS/Shendure_Lab/Projects/Many_targets_ideas/CROPt/191126_targATAC_full_1/210405_Tree42_ATAC_Analysis/bin5000upstreamPlusGeneBody/AllPlates_upstr5000_plusWholeGeneBody_cell_by_tss.txt"
ATAC_bins_tbl = read.table(ATAC_bins_tbl_name, sep = "\t", stringsAsFactors = FALSE)

##test total cell numbers
#reads_per_cell = rowSums(ATAC_bins_tbl)
#hist(log(reads_per_cell, base = 10), breaks = 100)
#read_cutoff = 0
#ATAC_bins_tbl = ATAC_bins_tbl[which(reads_per_cell > read_cutoff),]

#Add group info
ATAC_lin_group_file_name = "/Users/anna/Documents/UW_GS/Shendure_Lab/Projects/Many_targets_ideas/CROPt/190722_Tree29_other7plates/TARG_plates2thru8/200623_NewPipelineAnalysis/210405_ModTree42_ATAC_Analysis/All_Cell_ATAC_Newest_Groups"
ATAC_lin_groups = read.table(ATAC_lin_group_file_name, sep = "\t", stringsAsFactors = FALSE)
#ATAC_bins_tbl = ATAC_bins_tbl[which(rownames(ATAC_bins_tbl) %in% rownames(ATAC_lin_groups)),]
ATAC_lin_groups_ord = ATAC_lin_groups[match(rownames(ATAC_bins_tbl),rownames(ATAC_lin_groups)),]

##Add list of genes??
if (only_keep_set_list_of_genes == TRUE){
  tbl_w_list_of_genes_name = "/Users/anna/Documents/UW_GS/Shendure_Lab/Projects/Many_targets_ideas/CROPt/190722_Tree29_other7plates/TARG_plates2thru8/200623_NewPipelineAnalysis/210326_ModTree42_RNA_Analysis/ModTree42_HeatMaps/Matrix_File_W_Means"
  tbl_w_list_of_genes = read.table(tbl_w_list_of_genes_name, sep = "\t", stringsAsFactors = FALSE)
} 

if (just_gene_name_no_bins == TRUE){
  gene_metadata_with_names = read.table("/Users/anna/Documents/UW_GS/Shendure_Lab/Projects/Many_targets_ideas/CROPt/CROPt_pipeline_files/Commonly_sourced_files/unique_gencode.v29.gene_ids_gene_names_PLUS.bed", header = TRUE, stringsAsFactors = FALSE)
  gene_metadata_with_names$V2 = gene_metadata_with_names$undated_start_pos
  all_genes_on_chr = gene_metadata_with_names[,c(1,2,3,7)]
  
  ### SUBSET BY GENES IN FILE!!!!!!
  ATAC_bins_tbl = ATAC_bins_tbl[,which(colnames(ATAC_bins_tbl) %in% all_genes_on_chr$V7)]
}

ATAC_bins_tbl$Temp_col = ATAC_lin_groups_ord$V2

temp_for_counts = ATAC_bins_tbl[,c(1:(ncol(ATAC_bins_tbl)-1))]
counts_per_cell = rowSums(temp_for_counts)

cell_names_and_total_counts = as.data.frame(cbind(rownames(temp_for_counts), counts_per_cell))
colnames(cell_names_and_total_counts) = c("cell_name", "counts_per_cell")
write.table(cell_names_and_total_counts, "/Users/anna/Documents/UW_GS/Shendure_Lab/Projects/Many_targets_ideas/CROPt/190722_Tree29_other7plates/TARG_plates2thru8/200623_NewPipelineAnalysis/210405_ModTree42_ATAC_Analysis/211201_permutation_redo/ATAC_counts_per_cell", sep = "\t", row.names = FALSE, quote = FALSE)



##make a table with sums per bin per group
unique_groups = unique(ATAC_bins_tbl$Temp_col)
temp_ATAC_bins_tbl = ATAC_bins_tbl[,c(1:(ncol(ATAC_bins_tbl)-1))]
ATAC_bins_tbl_per_group = data.frame()
for(i in 1:length(unique_groups)){
  print(i)
  group = unique_groups[i]
  sub_tbl = ATAC_bins_tbl[which(ATAC_bins_tbl$Temp_col == group),c(1:(ncol(ATAC_bins_tbl)-1))]
  temp_bin_sums = colSums(sub_tbl)
  ATAC_bins_tbl_per_group = rbind(ATAC_bins_tbl_per_group, temp_bin_sums)
}
rownames(ATAC_bins_tbl_per_group) =paste0(rep("Group_", length(unique_groups)), unique_groups)

#ATAC_bins_tbl_per_group$group = unique_groups

ATAC_bins_tbl_per_group[1:10,1:10]


##NEW PSEUDOCOUNT TEST
#new_pseudocount = 1
#ATAC_bins_tbl_per_group = ATAC_bins_tbl_per_group + new_pseudocount
##NEW PSEUDOCOUNT TEST

##scale to group total read count 
total_read_count = rowSums(ATAC_bins_tbl_per_group)
hist(log(total_read_count, base = 10), breaks = 30)

scale_factor = median(total_read_count)
scale_factor

ATAC_bins_tbl_per_group_scaled = ATAC_bins_tbl_per_group/total_read_count*scale_factor
colnames(ATAC_bins_tbl_per_group_scaled) = colnames(ATAC_bins_tbl[,c(1:(ncol(ATAC_bins_tbl)-1))])

ATAC_bins_tbl_per_group_scaled[1:10,1:10]

##calculate means per bin
bins_n_means = as.data.frame(cbind(colMeans(ATAC_bins_tbl_per_group_scaled)))
hist(log(bins_n_means$V1, base = 2), breaks = 1000)
hist(log(bins_n_means$V1, base = 2), breaks = 100, xlim = range(0,4))

hist(bins_n_means$V1, breaks = 1000)

hist(bins_n_means$V1, breaks = 10000, xlim = range(0,50))

##subset to just higher bins
#min_mean_bin_ctoff = 70
min_mean_bin_ctoff = 5
bins_n_means_to_include =  as.data.frame(bins_n_means[which(bins_n_means$V1 > min_mean_bin_ctoff),,drop = FALSE])
ATAC_bins_tbl_per_group_scaled_high_bins = ATAC_bins_tbl_per_group_scaled[,which(colnames(ATAC_bins_tbl_per_group_scaled) %in% rownames(bins_n_means_to_include))]

ATAC_bins_tbl_per_group_scaled_high_bins = as.data.frame(t(t(ATAC_bins_tbl_per_group_scaled_high_bins)/as.numeric(bins_n_means_to_include$V1)))
hist(as.numeric(as.list(as.matrix(ATAC_bins_tbl_per_group_scaled_high_bins))), breaks = 100)

ATAC_bins_tbl_per_group_scaled_high_bins[1:10,1:10]

## add pseudocount and log
#pseudocount = 1 
pseudocount = 1 #NOT USED -- actually, maybe used
log_base = 2

##CHANGED THIS LINE FOR PSEUDOCOUNT TEST
log_ATAC_bins_tbl_per_group_scaled_high_bins = log(ATAC_bins_tbl_per_group_scaled_high_bins+pseudocount,base = log_base)
#log_ATAC_bins_tbl_per_group_scaled_high_bins = log(ATAC_bins_tbl_per_group_scaled_high_bins,base = log_base) 

log_ATAC_bins_tbl_per_group_scaled_high_bins$temp_groups = apply(as.data.frame(rownames(log_ATAC_bins_tbl_per_group_scaled_high_bins)), 1, split_on_delimiter, delim = "_", pos_to_return = 2)

log_ATAC_bins_tbl_per_group_scaled_high_bins$temp_groups = as.numeric(log_ATAC_bins_tbl_per_group_scaled_high_bins$temp_groups)

log_ATAC_bins_tbl_per_group_scaled_high_bins$temp_groups = as.numeric(log_ATAC_bins_tbl_per_group_scaled_high_bins$temp_groups)
log_ATAC_bins_tbl_per_group_scaled_high_bins = log_ATAC_bins_tbl_per_group_scaled_high_bins[order(log_ATAC_bins_tbl_per_group_scaled_high_bins$temp_groups),]
log_ATAC_bins_tbl_per_group_scaled_high_bins = log_ATAC_bins_tbl_per_group_scaled_high_bins[,c(1:(ncol(log_ATAC_bins_tbl_per_group_scaled_high_bins)-1))]

bin_logged_means = colMeans(log_ATAC_bins_tbl_per_group_scaled_high_bins)

log_ATAC_bins_tbl_per_group_scaled_high_bins = as.data.frame(t(t(log_ATAC_bins_tbl_per_group_scaled_high_bins)-bin_logged_means))
#hist(as.numeric(as.list(as.matrix(log_ATAC_bins_tbl_per_group_scaled_high_bins))), breaks = 100)

write.table(log_ATAC_bins_tbl_per_group_scaled_high_bins, paste0(Output_dir, "Matrix_File_W_Means"), sep = "\t", quote = FALSE)

add_chr_info = function(temp_mat_row, four_col_gene_file){
  temp_gene = as.character(temp_mat_row[1])
  temp_chr = as.character(four_col_gene_file[which(four_col_gene_file[,4] == temp_gene),1])
  return (temp_chr)
}

add_bin_info = function(temp_mat_row, four_col_gene_file){
  temp_gene = as.character(temp_mat_row[1])
  temp_bin = as.numeric(as.character(four_col_gene_file[which(four_col_gene_file[,4] == temp_gene),2]))
  return (temp_bin)
}

if(just_gene_name_no_bins == TRUE){
  gene_names_in_col_order = colnames(log_ATAC_bins_tbl_per_group_scaled_high_bins)
  df_gene_names_in_col_order = as.data.frame(gene_names_in_col_order)
  
  df_gene_names_in_col_order$chr = apply(as.data.frame(df_gene_names_in_col_order), 1, add_chr_info, four_col_gene_file = all_genes_on_chr)
  df_gene_names_in_col_order$bin = apply(as.data.frame(df_gene_names_in_col_order), 1, add_bin_info, four_col_gene_file = all_genes_on_chr)

  df_gene_names_in_col_order$concat = paste0(df_gene_names_in_col_order$gene_names_in_col_order, "_", df_gene_names_in_col_order$chr, "_", df_gene_names_in_col_order$bin)
  colnames(log_ATAC_bins_tbl_per_group_scaled_high_bins) = df_gene_names_in_col_order$concat
}

#Let's see if we can matrix this in a meaningful way...

temp_matrix = as.matrix(log_ATAC_bins_tbl_per_group_scaled_high_bins)
temp_matrix_melt = melt(temp_matrix)
temp_matrix_melt$width = 1
temp_matrix_melt$height = 1
hist(temp_matrix_melt$value, breaks = 100)

if(just_gene_name_no_bins == FALSE){
  temp_matrix_melt$chr = apply(as.data.frame(temp_matrix_melt)[,2, drop = FALSE], 1, split_on_delimiter, delim = "_", pos_to_return = 1)
  temp_matrix_melt$bin = apply(as.data.frame(temp_matrix_melt)[,2, drop = FALSE], 1, split_on_delimiter, delim = "_", pos_to_return = 2)
  
  if (encode_bins == FALSE){
    temp_matrix_melt$chr = gsub("X", "",temp_matrix_melt$chr)
  } else {
    temp_matrix_melt$chr = gsub("chr", "",temp_matrix_melt$chr)
    temp_matrix_melt[which(temp_matrix_melt$chr == "X"),]$chr = "23"
  }
  
  if (binned_by_gene == TRUE){
    temp_matrix_melt[which(temp_matrix_melt$chr == ""),]$chr = "23"
  }
  
  temp_matrix_melt$chr = as.numeric(as.character(temp_matrix_melt$chr))
  temp_matrix_melt = temp_matrix_melt[order(temp_matrix_melt$chr),]
}

if (just_gene_name_no_bins == TRUE){
  #add chr info
  temp_matrix_melt$chr = apply(as.data.frame(temp_matrix_melt)[,2, drop = FALSE], 1, split_on_delimiter, delim = "_", pos_to_return = 2)
  temp_matrix_melt$chr = gsub("chr", "",temp_matrix_melt$chr)
  temp_matrix_melt[which(temp_matrix_melt$chr == "X"),]$chr = "23"
  temp_matrix_melt$chr = as.numeric(temp_matrix_melt$chr)
  temp_matrix_melt = temp_matrix_melt[which(temp_matrix_melt$chr < 24),]
  #add start info as "bin"
  temp_matrix_melt$bin = apply(as.data.frame(temp_matrix_melt)[,2, drop = FALSE], 1, split_on_delimiter, delim = "_", pos_to_return = 3)
  temp_matrix_melt$chr = as.numeric(as.character(temp_matrix_melt$chr))
  temp_matrix_melt = temp_matrix_melt[order(as.numeric(as.character(temp_matrix_melt$bin))),]
  temp_matrix_melt = temp_matrix_melt[order(temp_matrix_melt$chr),]
}

temp_matrix_melt$group = apply(as.data.frame(temp_matrix_melt)[,1, drop = FALSE], 1, split_on_delimiter, delim = "_", pos_to_return = 2)
temp_matrix_melt$group = as.numeric(as.character(temp_matrix_melt$group))

temp_matrix_melt = temp_matrix_melt[order(temp_matrix_melt$group),]
hist(temp_matrix_melt$value, breaks = 1000, xlim = range(-2,2))

temp_matrix_melt_SAVED = temp_matrix_melt

temp_min = -.9
temp_max = .9

high_genes = "keep"
mean_or_median = "mean"

if (high_genes == "remove")
{
  ###THIS WILL NEED TO BE CHANGED!!
  temp_matrix_melt = temp_matrix_melt[which(!(as.character(temp_matrix_melt$Var1) %in% as.character(bad_genes_temp$Var1))),]
} else if (high_genes == "keep") {
  if(nrow(temp_matrix_melt[which(temp_matrix_melt$value > temp_max),]) > 0){
    temp_matrix_melt[which(temp_matrix_melt$value > temp_max),]$value = temp_max
  }
  if(nrow(temp_matrix_melt[which(temp_matrix_melt$value < temp_min),]) > 0){
    temp_matrix_melt[which(temp_matrix_melt$value < temp_min),]$value = temp_min
  }
}

temp_matrix_melt = as.data.frame(temp_matrix_melt)
temp_matrix_melt[which(temp_matrix_melt$chr == "23"),]$chr = "X"
temp_matrix_melt = temp_matrix_melt[which(temp_matrix_melt$chr != "24"),]

if (only_keep_set_list_of_genes == TRUE){
  temp_matrix_melt$gene = apply(as.data.frame(temp_matrix_melt)[,2, drop = FALSE], 1, split_on_delimiter, delim = "_", pos_to_return = 3)
  temp_matrix_melt = temp_matrix_melt[which(temp_matrix_melt$gene %in% rownames(tbl_w_list_of_genes)),]
}

if(just_gene_name_no_bins == TRUE){
  temp_matrix_melt$Var2 = apply(as.data.frame(temp_matrix_melt)[,2, drop = FALSE], 1, split_on_delimiter, delim = "_", pos_to_return = 1)
}

colnames(temp_matrix_melt) = c("Var2", "Var1", "value", "width", "height", "chr", "bin", "group")

colors_temp = rev(brewer.pal(11,"RdBu"))

List_of_plots = list()
plot_lengths_vect = vector()
chrom_names = unique(as.character(as.data.frame(temp_matrix_melt)$chr))

temp_matrix_melt = as.data.frame(temp_matrix_melt)
#temp_matrix_melt$Var1 = as.character(temp_matrix_melt$Var1)

for (i in 1:length(chrom_names))
{
  print (i)
  chrom = chrom_names[i]
  plot_title = paste0("chr",chrom)
  #gene_subset = all_genes_on_chr_in_data[which(all_genes_on_chr_in_data$V1 == chrom),]
  #temp_matrix_melt = as.data.frame(temp_matrix_melt)
  temp_matrix_melt_subset = temp_matrix_melt[which(as.character(temp_matrix_melt$chr) == chrom),]
  temp_matrix_melt_subset$Var1 = factor(temp_matrix_melt_subset$Var1, levels = unique(temp_matrix_melt_subset$Var1))
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
    theme_minimal() + ylim(rev(levels(temp_matrix_melt_subset$Var2))) +
    theme(legend.position="none") +
    theme(axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(title = plot_title)
  num_genes = length(unique(temp_matrix_melt_subset$Var1))
  pdf(paste0(Output_dir,paste0("Long_plot_NonLog_",mean_or_median,"_",temp_min,"to",temp_max,"_cutoff_",min_mean_bin_ctoff,"_",high_genes,"_", plot_title,".pdf"), sep = ""), width = (num_genes/5)+.5, height = 10)
  print(plot_individual)
  dev.off()
}
lay = t(matrix(plot_lengths_vect))
pdf(paste0(Output_dir,paste0("Long_plot_NonLog_",mean_or_median,"_",temp_min,"to",temp_max,"_cutoff_",min_mean_bin_ctoff,"_", high_genes,".pdf"), sep = ""), width = (150*(dim(lay)[2])/5000), height = (10))
#pdf(paste0(Output_dir,paste0("Opp_Long_plot_NonLog_",mean_or_median,"_",temp_min,"to",temp_max,"_cutoff_",mean_cutoff,".pdf"), sep = ""), width = 40, height = 10)
print(grid.arrange(grobs = List_of_plots, layout_matrix = lay, ncol = length(chrom_names)))
dev.off()




