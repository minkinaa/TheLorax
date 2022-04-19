#change lists_of_cells_vect, and maybe the main and output dirs, and run it.

library(ggplot2)
library(reshape2)
library(pheatmap)

generate_list_of_cells = function(list_of_ref_cells, lookup_table){
  temp_list_of_cells = vector()
  for(i in 1:length(list_of_ref_cells)){
    ref_cell = list_of_ref_cells[i]
    small_lookup = lookup_table[which(lookup_table$ref_cell == ref_cell),]
    temp_list_of_cells = c(temp_list_of_cells, small_lookup$cell_name)
  }
  return(temp_list_of_cells)
}

make_cells_split_into_lin_group_file = function(non_norm_mat_name, alg, upper_cutoff, PlusCounts){
  non_norm_mat = read.table(non_norm_mat_name, header = FALSE, sep = "\t")
  #print(head(non_norm_mat))
  rownames(non_norm_mat) = non_norm_mat$V1
  non_norm_mat = non_norm_mat[,c(2:ncol(non_norm_mat))]
  colnames(non_norm_mat) = rownames(non_norm_mat)
  sum_vect = apply(non_norm_mat, 1, sum)
  num_uniq = length(unique(sum_vect))
  #unique_rows_tbl = as.data.frame(table(non_norm_mat))
  print(paste0("Number_unique_rows: ", num_uniq))
  if(num_uniq > 1)
  {
    plot = pheatmap(non_norm_mat, clustering_method = alg, fontsize = 1, silent = TRUE)
    row_names_in_order = as.data.frame(rownames(non_norm_mat[plot$tree_row[["order"]],]))
    colnames(row_names_in_order) = c("cells_in_order")
    rownames(row_names_in_order) = row_names_in_order$cells_in_order
    
    cut_tree_1 = cutree(plot$tree_row, k = 1)
    cut_tree_df = as.data.frame(cut_tree_1)
    
    max_num_groups = nrow(non_norm_mat)
    for(i in 2:max_num_groups)
    {
      #col_name = paste0("cut_tree_", as.character(i))
      temp =  cutree(plot$tree_row, k = i)
      cut_tree_df = cbind(cut_tree_df, temp)
    }
    colnames(cut_tree_df) = c(1:max_num_groups)
    if (PlusCounts == TRUE)
    {
      write.table(cut_tree_df, file = paste0("Cells_PlusCounts_split_into_groups_", alg,"_AtOrAbove_", upper_cutoff,".txt"), quote = FALSE, sep = "\t")
    } else {
      write.table(cut_tree_df, file = paste0("Cells_split_into_groups_", alg,"_AtOrAbove_", upper_cutoff,".txt"), quote = FALSE, sep = "\t")
    }
  } else {
    cut_tree_df = non_norm_mat
    cut_tree_df$temp = 1 
    cut_tree_df$temp2 = 1 
    cut_tree_df = cut_tree_df[,c((ncol(cut_tree_df) - 1),ncol(cut_tree_df))]
    colnames(cut_tree_df) = c("1", "2")
    write.table(cut_tree_df, file = paste0("Cells_split_into_groups_", alg,"_AtOrAbove_", upper_cutoff,".txt"), quote = FALSE, sep = "\t")
  }
}

make_ref_cell_list = function(direct, vector_of_list_names){
  temp_list_of_cells = vector()
  for (i in 1:length(vector_of_list_names))
  {
    temp_cell_list_file = as.data.frame(read.table(paste0(direct,"/",vector_of_list_names[i]), sep = "\t", header = FALSE, stringsAsFactors = FALSE), drop = FALSE)
    print(head(temp_cell_list_file))
    temp_list = temp_cell_list_file$V1
    temp_list_of_cells = c(temp_list_of_cells, temp_list)
  }
  return(temp_list_of_cells)
}

corrected_cigar = TRUE #if change to False, may rewrite the TRUE files? Have not checked (keep this to TRUE)
recluster_cells = TRUE
from_list_of_cells = FALSE #will always be false for this script
curr_cutoff = 2 #change this -- this is the cutoff level to consider from 200219_make_LG_group_plots_for_combined_cell_groups.R 
alg = "ward.D2"
max_lineage_group_to_plot = 10 #script currently assumes we're workign w/ X10_LG___ : you can also make accommodations for other X's 

main_run_dir = "/Users/anna/Documents/UW_GS/Shendure_Lab/Projects/Many_targets_ideas/CROPt/190722_Tree29_other7plates/TARG_plates2thru8/200623_NewPipelineAnalysis/Corrected_CombinedPlates_match-5_uned-1/Combine_Like_Cells_200811_Take2"
colors_file = read.table("/Users/anna/Documents/UW_GS/Shendure_Lab/Projects/Many_targets_ideas/CROPt/190722_Tree29_other7plates/TARG_plates2thru8/200623_NewPipelineAnalysis/Corrected_CombinedPlates_match-5_uned-1/COLORS_Corrected_All_Plates_AMBcorr_Xcorr_cigar_match-5_uned-1", sep = "\t", comment.char = "", stringsAsFactors = FALSE)

cigar_files_dir = "/Users/anna/Documents/UW_GS/Shendure_Lab/Projects/Many_targets_ideas/CROPt/190722_Tree29_other7plates/TARG_plates2thru8/200623_NewPipelineAnalysis/Corrected_CombinedPlates_match-5_uned-1/"
Cigar_file = read.table(paste0(cigar_files_dir, "Corrected_All_Plates_AMBcorr_Xcorr_cigar_match-5_uned-1"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
Uncorr_Cigar_file = read.table(paste0(cigar_files_dir, "Uncorrected_All_Plates_AMBcorr_Xcorr_cigar_match-5_uned-1"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
Uncorr_Split_Target_Cigar_file = read.table(paste0(cigar_files_dir, "Uncorr_SplitTarget_All_Plates_AMBcorr_Xcorr_cigar_match-5_uned-1"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)

for (lg in 1:max_lineage_group_to_plot)
{
  curr_LG = lg
  print("Current Lingeage Group:", curr_LG)
  setwd(main_run_dir)
  dir_of_lists = paste0("Lineage_Plots_NoCounts_AtorAbove_",curr_cutoff)
  lists_of_cells_vect = c(paste0("CellList_X10_LG",curr_LG,"_cutoff",curr_cutoff,"_",alg))
  optional_cell_info_for_output = paste(lists_of_cells_vect, sep = "-")
  Output_dir = paste0("./", dir_of_lists, "/", paste("Gray_Partial_plots", optional_cell_info_for_output, sep = "-"), "/")
  system2("mkdir", args = Output_dir)
  
  new_list_of_cells = make_ref_cell_list(dir_of_lists, lists_of_cells_vect)
  
  for (ref_cell_num in 1:length(new_list_of_cells))
  {
    print(paste0("On cell number:", ref_cell_num))
    setwd(main_run_dir)
    #Output_dir = "/Users/anna/Documents/UW_GS/Shendure_Lab/Projects/Many_targets_ideas/CROPt/190722_Tree29_other7plates/TARG_plates2thru8/analysisNewPipeline/Fixed_data_all_plates/All_plates_combined/200219_combine_like_cells/Lineage_Plots_PlusCounts_AtorAbove_10/Partial_plots_p1A7_AATCATACGG/"
    if (from_list_of_cells == FALSE)
    {
      ref_cell_list = (new_list_of_cells[ref_cell_num])
      optional_cell_info = paste(ref_cell_list, sep = "-")
      #Output_dir = paste0("./", dir_of_lists, "/", paste("Partial_plots", optional_cell_info, sep = "-"), "/")
      #print(Output_dir)
      #system2("mkdir", args = Output_dir)
    }
    
    cell_lookup_tbl = read.table(paste0("Cell_to_RefCell_Lookup_Table_AtofAbove",curr_cutoff), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    #ref_cell_list = read.table(, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    #ref_cell_list = ("p1A7_AATCATACGG")
    if (corrected_cigar == TRUE)
    {
      Cigar_file = Cigar_file
      Uncorr_Cigar_file = Uncorr_Cigar_file
      Uncorr_Split_Target_Cigar_file = Uncorr_Split_Target_Cigar_file
      suffix = paste0("Corrected_",optional_cell_info)
    } else {
      Cigar_file = read.table("/Users/anna/Documents/UW_GS/Shendure_Lab/Projects/Many_targets_ideas/CROPt/190722_Tree29_other7plates/TARG_plates2thru8/analysisNewPipeline/2020_AllPlates/AllPlates_linRNA_Cigar_gr_0_ORDEREDwCOLNAMES", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
      suffix = paste0("Uncorrected_",optional_cell_info)
    }
    
    list_of_cells_to_plot = generate_list_of_cells(ref_cell_list, cell_lookup_tbl)
    partial_cigar = Cigar_file[which(rownames(Cigar_file) %in% list_of_cells_to_plot),]
    
    if (recluster_cells == FALSE)
    {
      ordered_cigar = partial_cigar
      ordered_cut_tree_df = ordered_cigar
      ordered_cut_tree_df$X10 = 1
      ordered_cut_tree_df = ordered_cut_tree_df[,c(ncol(ordered_cut_tree_df), ncol(ordered_cut_tree_df))]
    } else {
      setwd(Output_dir)
      file_name = paste0(paste0("Cigar_for_reclustering_NoColnames_", suffix))
      write.table(partial_cigar, file_name, col.names = FALSE, quote = FALSE, sep = "\t")
      program_name = "/Users/anna/Documents/UW_GS/Shendure_Lab/Projects/Many_targets_ideas/CROPt/CROPt_pipeline_files/Scripts/1900820_distance_matrix_bash_wrapper_script_2_run_num_as_argument.txt"
      #setwd(Output_dir)
      print(file_name)
      system2(program_name, args = c(file_name, as.character(curr_cutoff)))
      make_cells_split_into_lin_group_file(paste0("Non_Normalized_Distance_Mat_",curr_cutoff), alg, curr_cutoff, FALSE)
      cut_tree_df_first_set = read.table(paste0("Cells_split_into_groups_", alg,"_AtOrAbove_",curr_cutoff,".txt"), header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)
      for (i in (ncol(cut_tree_df_first_set)-1):2)
      {
        #print (i)
        cut_tree_df_first_set = cut_tree_df_first_set[order(cut_tree_df_first_set[,i]),]
      }
      ordered_cut_tree_df = cut_tree_df_first_set
      ordered_cigar = Cigar_file[match(rownames(ordered_cut_tree_df), rownames(Cigar_file)),]
      suffix = paste0(suffix,"_reclustered")
    }
    setwd(main_run_dir)
    
    if(ref_cell_num == 1)
    {
      merged_cigar = ordered_cigar
      merged_cigar_w_numbers = merged_cigar #add one w/ group_numbers
      merged_cigar_w_numbers$group_num = 1
      #print(paste0("NUM COLUMNS IN MERGED FILE 1:", ncol(merged_cigar_w_numbers)))
      group_counter = 2
      gray_vector = rep("gray",ncol(merged_cigar))
      merged_cigar = as.data.frame(rbind(ordered_cigar, gray_vector,gray_vector,gray_vector))
      print(nrow(merged_cigar))
      print(paste0("Length of merged cigar: ", nrow(merged_cigar)))
      
    } else {
      merged_cigar = as.data.frame(rbind(merged_cigar, ordered_cigar, gray_vector,gray_vector,gray_vector))
      ordered_cigar_w_num = ordered_cigar
      ordered_cigar_w_num$group_num = group_counter
      group_counter = group_counter+ 1
      #print(paste0("NUM COLUMNS IN MERGED FILE 2:", ncol(merged_cigar_w_numbers)))
      #print(paste0("NUM COLUMNS IN ORDERED FILE 2:", ncol(ordered_cigar_w_num)))
      merged_cigar_w_numbers = as.data.frame(rbind(merged_cigar_w_numbers, ordered_cigar_w_num))
      gray_vector = rep("gray",ncol(merged_cigar))
      print(paste0("Length of merged cigar: ", nrow(merged_cigar)))
    }
  }
  
  #write.table(merged_cigar_w_numbers, paste0(Output_dir,"X", Num_LGs, "_LG",j,"_", optional_cell_info_for_output,"UNCORRECTED.pdf", sep = ""))
  
  ordered_cigar = merged_cigar
  lineage_group_cutd_df = ordered_cigar
  lineage_group_cutd_df$temp1 = 1
  lineage_group_cutd_df$temp2 = 1
  lineage_group_cutd_df = lineage_group_cutd_df[,c(ncol(lineage_group_cutd_df)-1,ncol(lineage_group_cutd_df))]
  ordered_cut_tree_df = lineage_group_cutd_df
  
  #colors
  # all_editing_patterns_file = as.data.frame(read.table("/Users/anna/Documents/UW_GS/Shendure_Lab/Projects/Many_targets_ideas/CROPt/190722_Tree29_other7plates/TARG_plates2thru8/analysis/all_plates/test_folder_for_tree_calculations/All_Unique_Editing_Patters_5", header = FALSE), drop = FALSE)
  # num_unique = nrow(all_editing_patterns_file)
  # #set.seed(444)
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
  # cols[which(names(cols) == "gray")] = "black"
  #cols_3 = c(cols, "gray")
  #names(cols_3) = c(names(cols), "gray")
  #cols = cols_3
  #names(cols) = names(cols_3)
  
  cols = colors_file$V2
  names(cols) = colors_file$V1
  cols["gray"] = "black"
  
  ### ADDING A THING HERE WHERE YOU REPLACE CIGAR STRINGS W/ UNCORRECTED ONES!!
  temp_uncorr_cigar = Uncorr_Cigar_file[which(rownames(Uncorr_Cigar_file) %in% rownames(ordered_cigar)),]
  #gray_vector = rep("gray",ncol(temp_uncorr_cigar))
  gray_subset_of_ordered_cigar = ordered_cigar[which(!(rownames(ordered_cigar) %in% rownames(temp_uncorr_cigar))),]
  gray_subset_of_ordered_cigar = gray_subset_of_ordered_cigar[,c(1:ncol(temp_uncorr_cigar))]
  temp_uncorr_cigar = as.data.frame(rbind(temp_uncorr_cigar, gray_subset_of_ordered_cigar))
  temp_uncorr_cigar_ordered = temp_uncorr_cigar[match(rownames(ordered_cigar), rownames(temp_uncorr_cigar)),]
  
  ### 200804: Adding uncorrected but split target
  temp_uncorr_target_split_cigar = Uncorr_Split_Target_Cigar_file[which(rownames(Uncorr_Split_Target_Cigar_file) %in% rownames(ordered_cigar)),]
  gray_subset_of_ordered_cigar = ordered_cigar[which(!(rownames(ordered_cigar) %in% rownames(temp_uncorr_target_split_cigar))),]
  gray_subset_of_ordered_cigar = gray_subset_of_ordered_cigar[,c(1:ncol(temp_uncorr_target_split_cigar))]
  temp_uncorr_target_split_cigar = as.data.frame(rbind(temp_uncorr_target_split_cigar, gray_subset_of_ordered_cigar))
  temp_uncorr_target_split_cigar_ordered = temp_uncorr_target_split_cigar[match(rownames(ordered_cigar), rownames(temp_uncorr_target_split_cigar)),]
  
  #Make figure: 
  for (i in c(1))
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
        
        pdf(paste0(Output_dir,"X", Num_LGs, "_LG",j,"_", optional_cell_info_for_output,".pdf", sep = ""), width = 10, height = nrow(cigar_string_file_ordered)*.07 + 1.3)
        print(plot_tbl_ordered)
        dev.off()
        
        group_cells = as.data.frame(rownames(cigar_string_file_ordered), drop = FALSE)
        write.table(group_cells, file = paste0(Output_dir, "CellList_X", Num_LGs, "_LG",j, "_", optional_cell_info_for_output, sep = ""), row.names = FALSE, col.names = FALSE, quote = FALSE)
      
        write.table(merged_cigar_w_numbers, file = paste0(Output_dir, "CIGAR_W_GROUP_INFO_X", Num_LGs, "_LG",j, "_", optional_cell_info_for_output, sep = ""), row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")
        just_groups = c(rownames(merged_cigar_w_numbers), merged_cigar_w_numbers$group_num)
        
        write.table(just_groups, file = paste0(Output_dir, "LIST_OF_CELLS_W_GROUP_INFO_X", Num_LGs, "_LG",j, "_", optional_cell_info_for_output, sep = ""), row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
        
        ##ADDING IN A PART HERE TO PLOT THE UNEDITED STUFF; MAY BREAK EVERYTHIGN IN WHICH CASE DELETE EVERYTHING BELOW HERE.
        subset_cigar = temp_uncorr_cigar_ordered[which(rownames(temp_uncorr_cigar_ordered) %in% rownames(lineage_group_cutd_df)),]
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
        
        pdf(paste0(Output_dir,"X", Num_LGs, "_LG",j,"_", optional_cell_info_for_output,"UNCORRECTED.pdf", sep = ""), width = 10, height = nrow(cigar_string_file_ordered)*.07 + 1.3)
        print(plot_tbl_ordered)
        dev.off()
        
        ##  NOW ALSO ADDING A PART TO PLOT MAKE A SPLIT TARGET PLOT
        subset_cigar = temp_uncorr_target_split_cigar_ordered[which(rownames(temp_uncorr_target_split_cigar_ordered) %in% rownames(lineage_group_cutd_df)),]
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
        
        pdf(paste0(Output_dir,"X", Num_LGs, "_LG",j,"_", optional_cell_info_for_output,"UNCORRECTED_SPLIT_TARGET.pdf", sep = ""), width = 10, height = nrow(cigar_string_file_ordered)*.07 + 1.3)
        print(plot_tbl_ordered)
        dev.off()
      }
    }
  }
}






