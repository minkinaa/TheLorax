library(dplyr)
#library(Seurat)
library(ggplot2)
library(ggrepel)
library(reshape2)
library(robustbase) #colMedians
library("RColorBrewer")
library(gridExtra)
library(matrixStats)

calculate_num_more_extreme = function(row_of_numbers){
  temp_real_quotient = row_of_numbers[1]
  temp_ordered = row_of_numbers[order(row_of_numbers)]
  temp_pos = which(temp_ordered == temp_real_quotient)
  if(length(temp_pos) > 1){
    temp_middle = ceiling(length(temp_pos)/2)
    temp_pos = temp_pos[temp_middle]
  }
  return (temp_pos)
}

#Output_dir = "/Users/anna/Documents/UW_GS/Shendure_Lab/Projects/Many_targets_ideas/CROPt/190722_Tree29_other7plates/TARG_plates2thru8/200623_NewPipelineAnalysis/201117_All_RNA_Analysis/201125_Permuations_First_Attempt/"

args = commandArgs(trailingOnly=TRUE)
Output_dir = args[1]
gene_metadata_with_names_file = args[2]
group_file_name = args[3]
NonLog_Nor_counts_name = args[4]
Num_iterations = as.integer(args[5])
Run_num = as.integer(args[6])
cut_df_file = args[7]

#cut_df_file = "/Users/anna/Documents/UW_GS/Shendure_Lab/Projects/Many_targets_ideas/CROPt/190722_Tree29_other7plates/TARG_plates2thru8/200623_NewPipelineAnalysis/201117_All_RNA_Analysis/201120_Merge_Groups_First_Attempt/Consensus_Plot_New_Groups_Indicated/ReorderedCutDF_aka_cutdf"
cut_df = read.table(cut_df_file, sep = "\t", stringsAsFactors = FALSE)

#gene_metadata_with_names = read.table("/Users/anna/Documents/UW_GS/Shendure_Lab/Projects/Many_targets_ideas/CROPt/CROPt_pipeline_files/Commonly_sourced_files/unique_gencode.v29.gene_ids_gene_names_PLUS.bed", header = TRUE)
gene_metadata_with_names = read.table(gene_metadata_with_names_file, header = TRUE)
gene_metadata_with_names$V2 = gene_metadata_with_names$undated_start_pos
#chr_pos = read.table("/Users/anna/Documents/UW_GS/Shendure_Lab/Projects/Many_targets_ideas/CROPt/CROPt_pipeline_files/Commonly_sourced_files/chrom_length_midpoint_etc", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

#group_file = as.data.frame(read.table("/Users/anna/Documents/UW_GS/Shendure_Lab/Projects/Many_targets_ideas/CROPt/190722_Tree29_other7plates/TARG_plates2thru8/200623_NewPipelineAnalysis/201116_PostTreeBuild_AddAll_Cells/All_Cells_W_LinGroups_And_Rep_Cells_AND_NEW_Rep_Cells", header = TRUE, sep = "\t"), drop = FALSE)
group_file =  as.data.frame(read.table(group_file_name, sep = "\t", stringsAsFactors = FALSE))
#reorder = TRUE
#high_genes = "keep" # "remove" is other option

new_order = c(1:max(group_file$V2))
new_order_char = paste0("grp", new_order)

print("Reading in counts file...")

NonLog_Nor_counts = read.table(NonLog_Nor_counts_name, sep = "\t", row.names = 1)
NonLog_Nor_counts = NonLog_Nor_counts[match(rownames(group_file), rownames(NonLog_Nor_counts)), which(colnames(NonLog_Nor_counts) %in% as.character(gene_metadata_with_names$V7))]
###The above code gets rid of any cells which got tossed!

print("Done reading in counts file!")

NonLog_Nor_counts$Temp_col = group_file$V2

all_genes_on_chr = gene_metadata_with_names[,c(1,2,3,7)]
all_genes_on_chr_in_data = all_genes_on_chr[which(all_genes_on_chr$V7 %in% colnames(NonLog_Nor_counts)),]
#all_genes_on_chr_in_data = all_genes_on_chr_in_data[match(as.character(all_genes_on_chr_in_data$V7), colnames(NonLog_Nor_counts[,1:ncol(NonLog_Nor_counts)-1])),]
all_genes_on_chr_in_data = all_genes_on_chr_in_data[match(colnames(NonLog_Nor_counts[,1:ncol(NonLog_Nor_counts)-1]), as.character(all_genes_on_chr_in_data$V7)),]
all_genes_on_chr_in_data$mean_exp_in_data = colMeans(NonLog_Nor_counts[,1:ncol(NonLog_Nor_counts)-1])

pseudocount = 0.00000001

group_size_lookup_table = data.frame()

all_genes_plus_group_means = all_genes_on_chr_in_data[,1:4]
num_genes = nrow(all_genes_plus_group_means)

all_full_4_col_w_quotients = data.frame()

total_num_comparisons = 0

NonLog_Nor_counts = as.matrix(NonLog_Nor_counts)
start_time = proc.time()

time_since_last_iter = proc.time()

for(j in 2:ncol(cut_df)){
  print(paste0("Column #",j))
  unique_vals_in_previous_column = unique(cut_df[,j-1])
  for(k in 1:length(unique_vals_in_previous_column)){
    temp_group_num = unique_vals_in_previous_column[k]
    temp_cut_df_subset = cut_df[which(cut_df[,j-1] == temp_group_num),]
    list_of_unique_groups = unique(temp_cut_df_subset[,j])
    num_unique_groups_in_next_col = length(list_of_unique_groups)
    if (num_unique_groups_in_next_col > 1){
      ##make list of rep_cells in each group
      list_of_rep_cell_lists = list()
      for(m in 1:num_unique_groups_in_next_col){
        temp_goup = list_of_unique_groups[m]
        temp_list_of_rep_cells = rownames(temp_cut_df_subset[which(temp_cut_df_subset[,j] == temp_goup),])
        list_of_rep_cell_lists = c(list_of_rep_cell_lists, list(temp_list_of_rep_cells))
      }
      ##do pairwise comparisons
      for(n in 1:(length(list_of_rep_cell_lists)-1)){
        #print(paste0("n=", n))
        first_set_of_rep_cells = list_of_rep_cell_lists[n][[1]]
        first_group_col_and_group_num = paste0(j,"_",list_of_unique_groups[n])
        for (p in (n+1):length(list_of_rep_cell_lists)){
          #print(paste0("p=", p))
          second_set_of_rep_cells = list_of_rep_cell_lists[p][[1]]
          second_group_col_and_group_num = paste0(j,"_",list_of_unique_groups[p])
          print(paste0(first_group_col_and_group_num, " , ", second_group_col_and_group_num))
          ###do the comparisons here!!
          total_num_comparisons = total_num_comparisons + 1
          print (paste0("Time since last comparison:", proc.time() - time_since_last_iter))
          time_since_last_iter = proc.time() - time_since_last_iter
          print(paste0("On Comparison: ",total_num_comparisons))
          if(1 == 1){
            group_file_subset = group_file[which(group_file$new_rep_cell %in% c(first_set_of_rep_cells,second_set_of_rep_cells)),]
            group_file_subset$temp_group_num = 0
            group_file_subset[which(group_file_subset$new_rep_cell %in% first_set_of_rep_cells),]$temp_group_num = 1
            group_file_subset[which(group_file_subset$new_rep_cell %in% second_set_of_rep_cells),]$temp_group_num = 2
            
            num_cells_group_1 = nrow(group_file_subset[which(group_file_subset$temp_group_num == 1),])
            num_cells_group_2 = nrow(group_file_subset[which(group_file_subset$temp_group_num == 2),])
            temp_vect = as.data.frame(t(c(first_group_col_and_group_num, second_group_col_and_group_num, num_cells_group_1,num_cells_group_2)))
            group_size_lookup_table = rbind(group_size_lookup_table, temp_vect)
            
            temp_NonLog_Nor_counts = NonLog_Nor_counts[which(rownames(NonLog_Nor_counts) %in% rownames(group_file_subset)),]
            #temp_NonLog_Nor_counts = temp_NonLog_Nor_counts[match(rownames(group_file_subset), rownames(temp_NonLog_Nor_counts)),]
            temp_NonLog_Nor_counts = cbind(temp_NonLog_Nor_counts[,1:num_genes], group_file_subset$temp_group_num)
            #temp_NonLog_Nor_counts$Temp_col = group_file_subset$temp_group_num
            
            full_4_col_w_quotients = data.frame()
            
            for(iter in 1:Num_iterations){
              print(paste0("Comparison: ",total_num_comparisons, ", Iteration=",iter))
              #temp_4_col_tbl = data.frame()
              if(iter > 1){
                test_sample = sample(group_file_subset$temp_group_num,size = nrow(group_file_subset))
                temp_NonLog_Nor_counts = cbind(temp_NonLog_Nor_counts[,1:num_genes], test_sample)
                #temp_NonLog_Nor_counts$Temp_col = test_sample
              }
              temp_group_1_means = colMeans(temp_NonLog_Nor_counts[which(temp_NonLog_Nor_counts[,ncol(temp_NonLog_Nor_counts)] == 1),1:(ncol(temp_NonLog_Nor_counts)-1)])
              #print(paste0("Number of Means in grp1:", length(temp_group_1_means)))
              temp_group_2_means = colMeans(temp_NonLog_Nor_counts[which(temp_NonLog_Nor_counts[,ncol(temp_NonLog_Nor_counts)] == 2),1:(ncol(temp_NonLog_Nor_counts)-1)])
              #print(paste0("Number of Means in grp2:", length(temp_group_2_means)))
              temp_quotient = (temp_group_1_means+pseudocount)/(temp_group_2_means+pseudocount)
              #print(paste0("Number of Quotients:", length(temp_quotient)))
      
              if(iter == 1){
                first_set_sum_not_zero = colSums(temp_NonLog_Nor_counts[which(temp_NonLog_Nor_counts[,ncol(temp_NonLog_Nor_counts)] == 1),1:(ncol(temp_NonLog_Nor_counts)-1)] != 0)/nrow(temp_NonLog_Nor_counts[which(temp_NonLog_Nor_counts[,ncol(temp_NonLog_Nor_counts)] == 1),1:(ncol(temp_NonLog_Nor_counts)-1)])
                #print(paste0("Number of Sums in grp1:", length(first_set_sum_not_zero)))
                second_set_sum_not_zero = colSums(temp_NonLog_Nor_counts[which(temp_NonLog_Nor_counts[,ncol(temp_NonLog_Nor_counts)] == 2),1:(ncol(temp_NonLog_Nor_counts)-1)] != 0)/nrow(temp_NonLog_Nor_counts[which(temp_NonLog_Nor_counts[,ncol(temp_NonLog_Nor_counts)] == 2),1:(ncol(temp_NonLog_Nor_counts)-1)])
                #print(paste0("Number of Sums in grp2:", length(second_set_sum_not_zero)))
                #print(paste0("Num Genes: ", num_genes))
                
                full_4_col_w_quotients = as.data.frame(cbind(as.data.frame(rep(first_group_col_and_group_num, num_genes), stop = FALSE, stringsAsFactors = FALSE),
                                                             rep(second_group_col_and_group_num, num_genes, stop = FALSE, stringsAsFactors = FALSE), 
                                                             all_genes_plus_group_means$V7, temp_group_1_means, temp_group_2_means,
                                                             first_set_sum_not_zero, second_set_sum_not_zero,
                                                             temp_quotient), stringsAsFactors = FALSE, row.names = NULL)
                colnames(full_4_col_w_quotients) = c("Grp1_col_grpNum", "Grp2_col_grpNum", "Gene", "Grp1_Mean", "Grp2_Mean", "Grp1_Frac_NonZero","Grp2_Frac_NonZero", "Real_quotient")
              } else {
                full_4_col_w_quotients = cbind(full_4_col_w_quotients, temp_quotient)
              }
            }
            colnames(full_4_col_w_quotients) = c("Grp1_col_grpNum", "Grp2_col_grpNum", "Gene", "Grp1_Mean", "Grp2_Mean", "Grp1_Frac_NonZero","Grp2_Frac_NonZero", "Real_quotient", paste0("Perm_quotient", c(1:(Num_iterations-1))))
            full_4_col_w_quotients$pos_of_real = apply(full_4_col_w_quotients[,8:ncol(full_4_col_w_quotients)],1,calculate_num_more_extreme)
            all_full_4_col_w_quotients = rbind(all_full_4_col_w_quotients,full_4_col_w_quotients)
          }
        }
      }
    }
  }
}

end_time = proc.time()
print(end_time - start_time)

all_full_4_col_w_quotients = all_full_4_col_w_quotients[which(all_full_4_col_w_quotients$Grp1_Mean + all_full_4_col_w_quotients$Grp2_Mean != 0),]

rownames(all_full_4_col_w_quotients) = 1:nrow(all_full_4_col_w_quotients)

write.table(all_full_4_col_w_quotients, paste0(Output_dir, paste0("All_Pairwise_",Num_iterations,"_perm_Run", Run_num)), sep = "\t", quote = FALSE)
write.table(group_size_lookup_table, paste0(Output_dir, paste0("Group_Size_Lookup_",Num_iterations,"_perm_Run", Run_num)), sep = "\t", quote = FALSE)




