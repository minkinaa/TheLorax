
#### 200811 UPDATE: This is a modified version of 200225_combine_like_cells_for_loop 
#Noticed that rep cells were sometimes not representative of the group
#Would be good to potentially keep "rep cell" name to identify groups, but change the cigar string to be an average of groups

#this is a modified version of 200219_combine_like_cells w/ following mods:
#(1) for loop
#(2) internally runs scripts
#(3) has parametes for distances

#this is a great, super organized script which takes in a cigar file and reduces it down to non-unique cells
#there are places to run CPP scripts 
#you can keep increasing the upper_freq_cutoff from 2 to whatever. Really the whole thing should be a for loop which waits
#for you to run the cpp lines, but this totally works for now. 


library(pheatmap)

edits_to_string = function(table_row){
  #return("idk")
  temp_row = table_row
  temp_row = gsub("70M:","S", temp_row)
  temp_row = gsub("X","S", temp_row)
  temp_row = gsub("AMB","S", temp_row)
  temp_row_string = paste(temp_row, collapse = "")
  return(temp_row_string)
}

generate_cpp_call = function(map_to_cells, LGs_of_map_to_cells, map_cells, cutoff)
{
  print(paste0("/Users/anna/Documents/UW_GS/Shendure_Lab/Projects/Many_targets_ideas/CROPt/CROPt_pipeline_files/Scripts/191017_wrapper_script_for_dist_cpp_script_run_num_as_arg.txt ", map_to_cells," ", LGs_of_map_to_cells, " ", map_cells, " ", cutoff))
}

modify_edit_string = function(edit_cell_tbl, LinG_above_cutoff, LinG_below_cutoff, edit_cell_tbl_2){
  if(!(edit_cell_tbl[2] %in% rownames(LinG_below_cutoff))) #if not in unique edits
  {
    return(as.character(edit_cell_tbl[1]))
  }
  else
  {
    temp_cell_name = edit_cell_tbl[2] #name of current cell
    temp_lg = as.numeric(as.character(LinG_below_cutoff[which(rownames(LinG_below_cutoff) == temp_cell_name),1]))
    old_cell_name = as.character(rownames(LinG_above_cutoff[which(LinG_above_cutoff$LG == temp_lg),]))
    old_cell_edit_string = as.character(edit_cell_tbl_2[which(rownames(edit_cell_tbl_2) == old_cell_name),]$edit_string)
    return(old_cell_edit_string)
    #return(paste0("wait is this something? ",old_cell_name))
  }
}

find_cell_w_matching_string = function(all_cell_file, unique_cell_file)
{
  temp_concat_cigar = all_cell_file[2]
  temp_matching_cell_name = as.character(unique_cell_file[which(unique_cell_file$edit_string == temp_concat_cigar),1])
  return(temp_matching_cell_name)
}

generate_cpp_call_to_calc_dist_matrix_NoCounts = function(upper_freq_cutoff) {
  print(paste("/Users/anna/Documents/UW_GS/Shendure_Lab/Projects/Many_targets_ideas/CROPt/CROPt_pipeline_files/Scripts/1900820_distance_matrix_bash_wrapper_script_2_run_num_as_argument.txt", paste0("PostDistribution_Rep_Cell_Cigar_file_AtorAbove",upper_freq_cutoff, "_noColNames"), upper_freq_cutoff))
}

generate_cpp_call_to_calc_dist_matrix_PlusCounts = function(upper_freq_cutoff) {
  print(paste("/Users/anna/Documents/UW_GS/Shendure_Lab/Projects/Many_targets_ideas/CROPt/CROPt_pipeline_files/Scripts/1900820_distance_matrix_bash_wrapper_script_2_run_num_as_argument.txt", paste0("PostDistribution_Rep_Cell_Cigar_PlusCounts_AtorAbove",upper_freq_cutoff, "_noColNames"), paste0(upper_freq_cutoff,"_plusCounts")))
}

make_cells_split_into_lin_group_file = function(non_norm_mat_name, alg, upper_cutoff, PlusCounts){
  non_norm_mat = read.table(non_norm_mat_name, header = FALSE, sep = "\t")
  #print(head(non_norm_mat))
  rownames(non_norm_mat) = non_norm_mat$V1
  non_norm_mat = non_norm_mat[,c(2:ncol(non_norm_mat))]
  colnames(non_norm_mat) = rownames(non_norm_mat)
  
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
}

generate_new_folder_call = function(cutoff, PlusCounts) {
  if (PlusCounts == TRUE)
  {
    print(paste0("mkdir Lineage_Plots_PlusCounts_AtorAbove_",cutoff))
  } else {
    print(paste0("mkdir Lineage_Plots_NoCounts_AtorAbove_",cutoff))
  }
}



###   START HERE   ####

num_targets = 33

#setwd("/Users/anna/Documents/UW_GS/Shendure_Lab/Projects/Many_targets_ideas/CROPt/190722_Tree29_other7plates/TARG_plates2thru8/analysisNewPipeline/Fixed_data_all_plates/All_plates_combined_match5_uned1/200226_combine_like_cells")
setwd("/Users/anna/Documents/UW_GS/Shendure_Lab/Projects/Many_targets_ideas/CROPt/190722_Tree29_other7plates/TARG_plates2thru8/200623_NewPipelineAnalysis/Corrected_CombinedPlates_match-5_uned-1/Combine_Like_Cells_200811_Take2")
system2("touch", args = c("Log_File"))
system2("echo", args = c(">Log_File"))

#table contains AMB and X-corrected versions of all cells
all_cell_tbl = read.table("/Users/anna/Documents/UW_GS/Shendure_Lab/Projects/Many_targets_ideas/CROPt/190722_Tree29_other7plates/TARG_plates2thru8/200623_NewPipelineAnalysis/Corrected_CombinedPlates_match-5_uned-1/Corrected_All_Plates_AMBcorr_Xcorr_cigar_match-5_uned-1", sep = "\t", stringsAsFactors = FALSE, header = TRUE)
all_cell_tbl_mod = all_cell_tbl #the version of the table to add columns to 

#all_cell_tbl_mod$edit_string = apply(all_cell_tbl_mod, 1, edits_to_string) #make edit string


## COME BACK HERE WHEN RESTARTING THE LOOP!! 

highest_cutoff_value = 10
weights_same_edit = 5
weights_both_70 = 1
system2("echo", args = c(paste0("weights_same_edit: ",weights_same_edit), ">>", "Log_File"))
system2("echo", args = c(paste0("weights_both_70: ",weights_same_edit), ">>", "Log_File"))

for (i in 2:highest_cutoff_value){
  upper_freq_cutoff = i
  print(paste0("We are at round number: ", upper_freq_cutoff))
  if (upper_freq_cutoff == 2)
  {
    all_cell_tbl_mod$edit_string = apply(all_cell_tbl_mod, 1, edits_to_string) #make initial edit string
  } else {
    all_cell_tbl_mod$edit_string = edit_string_and_cell_name_tbl$new_edit_string
  }
  
  edit_string_freq_tbl = as.data.frame(table(all_cell_tbl_mod$edit_string))
  
  #make a table w/ a single cell representing each edit_string
  
  ##new stuff 200811. Here, we are only using "representative" cells whose actual cigar string matches rep. cigar string exactly.
  temp_all_cell_tbl_mod = all_cell_tbl_mod
  ##print ("is the problem here? 1")
  temp_all_cell_tbl_mod$old_cigar_string = apply(temp_all_cell_tbl_mod[,c(1:(ncol(temp_all_cell_tbl_mod)-1))], 1, edits_to_string)
  ##print ("is the problem here? 2")
  temp_all_cell_tbl_mod_sub = temp_all_cell_tbl_mod[which(temp_all_cell_tbl_mod$edit_string == temp_all_cell_tbl_mod$old_cigar_string),]
  #####
  rep_cell_cigar_tbl = temp_all_cell_tbl_mod_sub[match(unique(temp_all_cell_tbl_mod_sub$edit_string), temp_all_cell_tbl_mod_sub$edit_string),c(1:(ncol(temp_all_cell_tbl_mod_sub)-1))]
  ##line above was modified, too.
  
  edit_string_below_cutoff = edit_string_freq_tbl[which(edit_string_freq_tbl$Freq < upper_freq_cutoff),]
  edit_string_above_cutoff = edit_string_freq_tbl[which(edit_string_freq_tbl$Freq >= upper_freq_cutoff),]
  
  rep_cell_cigar_tbl_below_cutoff = all_cell_tbl_mod[which(all_cell_tbl_mod$edit_string %in% edit_string_below_cutoff$Var1),c(1:num_targets)]
  #rep_cell_cigar_tbl_below_cutoff = rep_cell_cigar_tbl[which(rep_cell_cigar_tbl$edit_string %in% edit_string_below_cutoff$Var1),c(1:num_targets)]
  rep_cell_cigar_tbl_above_cutoff = rep_cell_cigar_tbl[which(rep_cell_cigar_tbl$edit_string %in% edit_string_above_cutoff$Var1),c(1:num_targets)]
  
  lineage_group_file = rep_cell_cigar_tbl_above_cutoff
  lineage_group_file$LG = c(1:nrow(lineage_group_file))
  lineage_group_file = lineage_group_file[,c(ncol(lineage_group_file), ncol(lineage_group_file))]
  
  #setwd("/Users/anna/Documents/UW_GS/Shendure_Lab/Projects/Many_targets_ideas/CROPt/190722_Tree29_other7plates/TARG_plates2thru8/analysisNewPipeline/Fixed_data_all_plates/All_plates_combined/200219_combine_like_cells")
  write.table(rep_cell_cigar_tbl_below_cutoff, paste0("Cigar_Rep_Cells_Below_",upper_freq_cutoff), sep = "\t", col.names = TRUE, quote = FALSE)
  write.table(rep_cell_cigar_tbl_above_cutoff, paste0("Cigar_Rep_Cells_AtorAbove_",upper_freq_cutoff), sep = "\t", col.names = TRUE, quote = FALSE)
  
  write.table(lineage_group_file, paste0("LG_Rep_Cells_AtorAbove_",upper_freq_cutoff), sep = "\t", col.names = TRUE, quote = FALSE)
  
  #RUN CPP TO FIND NEAREST NEIGHBOR CELLS
  #generate_cpp_call(paste0("Cigar_Rep_Cells_AtorAbove_",upper_freq_cutoff), paste0("LG_Rep_Cells_AtorAbove_",upper_freq_cutoff), paste0("Cigar_Rep_Cells_Below_",upper_freq_cutoff), upper_freq_cutoff)
  cpp_prog_name = "/Users/anna/Documents/UW_GS/Shendure_Lab/Projects/Many_targets_ideas/CROPt/CROPt_pipeline_files/Scripts/191017_wrapper_script_for_dist_cpp_script_run_num_as_arg_Weights_As_Args.txt"
  system2(cpp_prog_name, args = c(paste0("Cigar_Rep_Cells_AtorAbove_",upper_freq_cutoff), paste0("LG_Rep_Cells_AtorAbove_",upper_freq_cutoff), paste0("Cigar_Rep_Cells_Below_",upper_freq_cutoff), upper_freq_cutoff, weights_same_edit, weights_both_70))
  
  #Use lineage group assignments to change unique edit_strings to common edit_strings
  new_lineage_groups = read.table(paste0("AllCells_Lineage_Groups_",upper_freq_cutoff))
  New_LG_above_cutoff = new_lineage_groups[which(rownames(new_lineage_groups) %in% rownames(rep_cell_cigar_tbl_above_cutoff)),]
  New_LG_below_cutoff = new_lineage_groups[which(rownames(new_lineage_groups) %in% rownames(rep_cell_cigar_tbl_below_cutoff)),]
  
  #correct cigar strings
  #we'd like to end up with a table the structure of all_cell_tbl_mod (w/ edit string)
  
  edit_string_and_cell_name_tbl = all_cell_tbl_mod
  edit_string_and_cell_name_tbl$cell_name = rownames(edit_string_and_cell_name_tbl)
  edit_string_and_cell_name_tbl = edit_string_and_cell_name_tbl[,c(ncol(edit_string_and_cell_name_tbl)-1, ncol(edit_string_and_cell_name_tbl))]
  
  edit_string_and_cell_name_tbl$new_edit_string = apply(edit_string_and_cell_name_tbl, 1, modify_edit_string, LinG_above_cutoff = New_LG_above_cutoff, LinG_below_cutoff = New_LG_below_cutoff, edit_cell_tbl_2 = edit_string_and_cell_name_tbl)
  #table(as.data.frame(table(edit_string_and_cell_name_tbl$new_edit_string))$Freq)
  
  #make a new table of representative cells for input into distance script (for plotting)
  cigar_string_tbl_w_new_edits = all_cell_tbl_mod
  cigar_string_tbl_w_new_edits$edit_string = edit_string_and_cell_name_tbl$new_edit_string
  
  ### 200811 -- OK THERE MIGHT BE AN ISSUE BELOW:
  ##adding new stuff here, too:
  
  temp_all_cell_tbl_mod = cigar_string_tbl_w_new_edits
  ##print ("is the problem here? 1")
  temp_all_cell_tbl_mod$old_cigar_string = apply(temp_all_cell_tbl_mod[,c(1:(ncol(temp_all_cell_tbl_mod)-1))], 1, edits_to_string)
  ##print ("is the problem here? 2")
  temp_all_cell_tbl_mod_sub = temp_all_cell_tbl_mod[which(temp_all_cell_tbl_mod$edit_string == temp_all_cell_tbl_mod$old_cigar_string),]
  #####
  #rep_cell_cigar_tbl = temp_all_cell_tbl_mod_sub[match(unique(temp_all_cell_tbl_mod_sub$edit_string), temp_all_cell_tbl_mod_sub$edit_string),c(1:(ncol(temp_all_cell_tbl_mod_sub)-1))]
  ###
  
  new_rep_cell_cigar_tbl = temp_all_cell_tbl_mod_sub[match(unique(temp_all_cell_tbl_mod_sub$edit_string), temp_all_cell_tbl_mod_sub$edit_string),c(1:(ncol(temp_all_cell_tbl_mod_sub)-1))]
  new_edit_strings_and_counts = as.data.frame(table(edit_string_and_cell_name_tbl$new_edit_string))
  new_edit_strings_and_counts = new_edit_strings_and_counts[match(new_rep_cell_cigar_tbl$edit_string, new_edit_strings_and_counts$Var1),]
  
  new_rep_cell_cigar_tbl$cell_count = new_edit_strings_and_counts$Freq
  new_rep_cell_cigar_tbl$cell_name_and_count = paste0(rownames(new_rep_cell_cigar_tbl),"-",as.character(new_rep_cell_cigar_tbl$cell_count))
  new_rep_cell_cigar_tbl_NoCounts = new_rep_cell_cigar_tbl[,c(1:num_targets)]
  new_rep_cell_cigar_tbl_PlusCounts = new_rep_cell_cigar_tbl
  rownames(new_rep_cell_cigar_tbl_PlusCounts) = new_rep_cell_cigar_tbl_PlusCounts$cell_name_and_count
  new_rep_cell_cigar_tbl_PlusCounts = new_rep_cell_cigar_tbl_PlusCounts[,c(1:num_targets)]
  
  write.table(new_rep_cell_cigar_tbl_NoCounts, paste0("PostDistribution_Rep_Cell_Cigar_file_AtorAbove",upper_freq_cutoff), sep = "\t", col.names = TRUE, quote = FALSE)
  write.table(new_rep_cell_cigar_tbl_NoCounts, paste0("PostDistribution_Rep_Cell_Cigar_file_AtorAbove",upper_freq_cutoff, "_noColNames"), sep = "\t", col.names = FALSE, quote = FALSE)
  
  write.table(new_rep_cell_cigar_tbl_PlusCounts, paste0("PostDistribution_Rep_Cell_Cigar_PlusCounts_AtorAbove",upper_freq_cutoff), sep = "\t", col.names = TRUE, quote = FALSE)
  write.table(new_rep_cell_cigar_tbl_PlusCounts, paste0("PostDistribution_Rep_Cell_Cigar_PlusCounts_AtorAbove",upper_freq_cutoff, "_noColNames"), sep = "\t", col.names = FALSE, quote = FALSE)
  
  #Make lookup table
  all_cell_name_edit = edit_string_and_cell_name_tbl[,c(ncol(edit_string_and_cell_name_tbl)-1, ncol(edit_string_and_cell_name_tbl))]
  unique_cell_name_edit = new_rep_cell_cigar_tbl
  unique_cell_name_edit$cell_name = rownames(unique_cell_name_edit)
  unique_cell_name_edit = unique_cell_name_edit[,c(ncol(unique_cell_name_edit), (ncol(unique_cell_name_edit) - 3))]
  
  all_cell_name_edit$ref_cell = apply(all_cell_name_edit, 1, find_cell_w_matching_string, unique_cell_file = unique_cell_name_edit)
  
  look_up_tbl = all_cell_name_edit[,c(1,3)]
  write.table(look_up_tbl, paste0("Cell_to_RefCell_Lookup_Table_AtofAbove", upper_freq_cutoff), row.names = FALSE, sep = "\t", quote = FALSE)
  
  #Ok let's make some files to make distance matrices for plotting 
  #RUN CPP FILE FROM CALLS GENERATED BELOW
  cpp_program_name_2 = "/Users/anna/Documents/UW_GS/Shendure_Lab/Projects/Many_targets_ideas/CROPt/CROPt_pipeline_files/Scripts/1900820_distance_matrix_bash_wrapper_script_2_run_num_as_argument.txt"
  system2(cpp_program_name_2, args = c(paste0("PostDistribution_Rep_Cell_Cigar_file_AtorAbove",upper_freq_cutoff,"_noColNames"), upper_freq_cutoff))
  system2(cpp_program_name_2, args = c(paste0("PostDistribution_Rep_Cell_Cigar_PlusCounts_AtorAbove",upper_freq_cutoff,"_noColNames"), paste0(upper_freq_cutoff, "_plusCounts")))
  #generate_cpp_call_to_calc_dist_matrix_NoCounts(upper_freq_cutoff)
  #generate_cpp_call_to_calc_dist_matrix_PlusCounts(upper_freq_cutoff)
  
  #Let's make some cells split into group files!
  make_cells_split_into_lin_group_file(paste0("Non_Normalized_Distance_Mat_",upper_freq_cutoff, "_plusCounts"), "ward.D2", upper_freq_cutoff, TRUE)
  make_cells_split_into_lin_group_file(paste0("Non_Normalized_Distance_Mat_",upper_freq_cutoff), "ward.D2", upper_freq_cutoff, FALSE)
  
  #bash calls to make folders for lineage group files 
  system2("mkdir", args = c(paste0("Lineage_Plots_PlusCounts_AtorAbove_", upper_freq_cutoff)))
  system2("mkdir", args = c(paste0("Lineage_Plots_NoCounts_AtorAbove_", upper_freq_cutoff)))
  #generate_new_folder_call(upper_freq_cutoff, TRUE)
  #generate_new_folder_call(upper_freq_cutoff, FALSE)
}
# RESTART THE LOOP!










