#This is dope
#change max cutoff, withCounts, and current directory and run! 
#This script follows directly after 200225_combine_like_cells_for_loop.R

library(ggplot2)
library(reshape2)

setwd("/Users/anna/Documents/UW_GS/Shendure_Lab/Projects/Many_targets_ideas/CROPt/190722_Tree29_other7plates/TARG_plates2thru8/200623_NewPipelineAnalysis/Corrected_CombinedPlates_match-5_uned-1/Combine_Like_Cells_200811_Take2")

max_cutoff = 10
withCounts = FALSE

colors_file = read.table("/Users/anna/Documents/UW_GS/Shendure_Lab/Projects/Many_targets_ideas/CROPt/190722_Tree29_other7plates/TARG_plates2thru8/200623_NewPipelineAnalysis/Corrected_CombinedPlates_match-5_uned-1/COLORS_Corrected_All_Plates_AMBcorr_Xcorr_cigar_match-5_uned-1", sep = "\t", comment.char = "", stringsAsFactors = FALSE)

for (i in 2:max_cutoff)
{
  curr_cutoff = i
  print(paste0("Currently on cutoff: ", curr_cutoff))
  #withCounts = FALSE
  alg = "ward.D2"
  suffix = paste0("cutoff",curr_cutoff,"_",alg)
  
  if(withCounts == TRUE){
    Output_dir = paste0("Lineage_Plots_PlusCounts_AtorAbove_",curr_cutoff,"/")
    Cigar_file = read.table(paste0("PostDistribution_Rep_Cell_Cigar_PlusCounts_AtorAbove", curr_cutoff), sep = "\t", row.names = 1, header = TRUE)
    cut_tree_df_first_set = read.table(paste0("Cells_PlusCounts_split_into_groups_", alg,"_AtOrAbove_",curr_cutoff,".txt"), header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)
  } else{
    Output_dir = paste0("Lineage_Plots_NoCounts_AtorAbove_",curr_cutoff,"/")
    Cigar_file = read.table(paste0("PostDistribution_Rep_Cell_Cigar_file_AtorAbove", curr_cutoff), sep = "\t", row.names = 1, header = TRUE)
    cut_tree_df_first_set = read.table(paste0("Cells_split_into_groups_", alg,"_AtOrAbove_",curr_cutoff,".txt"), header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)
  }
  
  for (i in (ncol(cut_tree_df_first_set)-1):2)
  {
    print (i)
    cut_tree_df_first_set = cut_tree_df_first_set[order(cut_tree_df_first_set[,i]),]
  }
  
  ordered_cut_tree_df = cut_tree_df_first_set
  
  ordered_cigar = Cigar_file[match(rownames(ordered_cut_tree_df), rownames(Cigar_file)),]
  
  # all_editing_patterns_file = as.data.frame(read.table("/Users/anna/Documents/UW_GS/Shendure_Lab/Projects/Many_targets_ideas/CROPt/190722_Tree29_other7plates/TARG_plates2thru8/analysis/all_plates/test_folder_for_tree_calculations/All_Unique_Editing_Patters_5", header = FALSE, stringsAsFactors = FALSE), drop = FALSE)
  # 
  # num_unique = nrow(all_editing_patterns_file)
  # set.seed(1)
  # cols = rainbow(num_unique, s=.6, v=.9)[sample(1:num_unique,num_unique)]
  # all_editing_patterns_file$colors = cols
  # 
  # all_editing_patterns_file[which(all_editing_patterns_file$V1 == "70M:"),]$colors = "#FFFFFF"
  # all_editing_patterns_file[which(all_editing_patterns_file$V1 == "X"),]$colors = "black"
  # all_editing_patterns_file[which(all_editing_patterns_file$V1 == "AMB"),]$colors = "red"
  # cols = as.character(all_editing_patterns_file$colors)
  # names(cols) = as.character(all_editing_patterns_file$V1)
  # 
  # unique_edits_in_new_set = as.data.frame(unique(as.vector(as.matrix(ordered_cigar))))
  # colnames(unique_edits_in_new_set) = "V1"
  # unique_edits_in_new_set$V1 = as.character(unique_edits_in_new_set$V1)
  # just_new = as.data.frame(unique_edits_in_new_set[which(!(unique_edits_in_new_set$V1 %in% all_editing_patterns_file$V1)),], drop = FALSE)
  # colnames(just_new) = "V1"
  # set.seed(2)
  # cols_2 = rainbow(num_unique+nrow(just_new), s=.6, v=.9)[sample(num_unique+nrow(just_new),num_unique+nrow(just_new))]
  # new_colors = cols_2[which(!(cols_2 %in% cols))]
  # 
  # cols = c(cols, cols_2[1:nrow(just_new)])
  # names(cols) = c(as.character(all_editing_patterns_file$V1),as.character(just_new$V1))
  # 
  
  cols = colors_file$V2
  names(cols) = colors_file$V1
  
  for (i in c(1,7,10,20))
  {
    Num_LGs = i
    for (j in 1:Num_LGs)
    {
      lineage_group_cutd_df = ordered_cut_tree_df[which(ordered_cut_tree_df[,i] == j),]
      print(nrow(lineage_group_cutd_df ))
      {
        subset_cigar = ordered_cigar[which(rownames(ordered_cigar) %in% rownames(lineage_group_cutd_df)),]
        #print(nrow(subset_cigar))
        cigar_string_file_ordered = subset_cigar
        cigar_string_file_ordered_mat = as.matrix(cigar_string_file_ordered)
        cigar_string_file_ordered_mat_melt = melt(cigar_string_file_ordered_mat)
        cigar_string_file_ordered_mat_melt$value = as.factor(cigar_string_file_ordered_mat_melt$value)
        cigar_string_file_ordered_mat_melt$width = 1
        cigar_string_file_ordered_mat_melt$height = .8
        #print(head(cigar_string_file_ordered_mat_melt))
        
        plot_tbl_ordered = ggplot(cigar_string_file_ordered_mat_melt, aes(x = Var2, y = Var1)) +
          geom_tile(aes(fill=value, width = cigar_string_file_ordered_mat_melt$width, height = cigar_string_file_ordered_mat_melt$height))
        
        plot_tbl_ordered = plot_tbl_ordered + scale_fill_manual(values = cols) + 
          theme_minimal() + theme(legend.position="none") + 
          theme(axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x=element_text(angle=90,hjust=1), axis.text.y=element_text(size = 6)) +
          labs(tag = j) +
          theme(plot.tag.position = c(.01, 4.83/(nrow(cigar_string_file_ordered) + 6))) +
          #theme(plot.tag.position = "bottomleft") +
          theme(plot.tag = element_text(size = 50, vjust = 0, hjust = 0))
        
        pdf(paste0(Output_dir,"X", Num_LGs, "_LG",j,"_", suffix,".pdf", sep = ""), width = 10, height = nrow(cigar_string_file_ordered)*.07 + 1.3)
        print(plot_tbl_ordered)
        dev.off()
        
        group_cells = as.data.frame(rownames(cigar_string_file_ordered), drop = FALSE)
        write.table(group_cells, file = paste0(Output_dir, "CellList_X", Num_LGs, "_LG",j, "_", suffix, sep = ""), row.names = FALSE, col.names = FALSE, quote = FALSE)
      }
    }
  }
}
