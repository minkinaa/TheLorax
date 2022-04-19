#In this script, we'll attempt to make a tree file out of the output from: 201106_building_a_tree_2.cpp
# output called target_edits_over_time_outfile in cpp script

args = commandArgs(trailingOnly=TRUE)
output_dir = args[1]
input_tree_file = args[2]

#output_dir = "/Users/anna/Documents/UW_GS/Shendure_Lab/Projects/Many_targets_ideas/CROPt/190722_Tree29_other7plates/TARG_plates2thru8/200623_NewPipelineAnalysis/210322_ModifyingTree/210323_modified_tree_files/"

#input_tree_file = "/Users/anna/Documents/UW_GS/Shendure_Lab/Projects/Many_targets_ideas/CROPt/190722_Tree29_other7plates/TARG_plates2thru8/200623_NewPipelineAnalysis/201028_TreeAttempt1/Corrected_p1_p3/With_Added_group_splits/201109_tree_splitting_second_attempt/target_edits_over_time_outfile"
input_tree = read.table(input_tree_file, sep = "\t", stringsAsFactors = FALSE)
colnames(input_tree) = as.character(input_tree[1,])
input_tree = input_tree[2:nrow(input_tree),]

#Add one to target pos so correspond to actual positions
input_tree$TargetPos = as.numeric(input_tree$TargetPos)
input_tree$TargetPos = input_tree$TargetPos + 1

input_tree$ChangeNum = as.numeric(input_tree$ChangeNum)

input_tree$ChangeNum_Target_Pos_Edit = paste0(input_tree$ChangeNum, "_", input_tree$TargetPos, "_", input_tree$Edit)

#Here, we'll make a vector of edits in order for each group
changes_df = as.data.frame(matrix(data = "NONE", nrow = length(unique(input_tree$Group)), ncol = max(input_tree$ChangeNum)), stringsAsFactors = FALSE)
rownames(changes_df) = unique(input_tree$Group)

for(i in 1:nrow(changes_df)){
  cell_name = rownames(changes_df[i,])
  sub_input_tree = input_tree[which(input_tree$Group == cell_name),]
  sub_input_tree = sub_input_tree[order(sub_input_tree$ChangeNum),]
  last_edit = "NONE"
  for (j in 1:nrow(sub_input_tree)){
    print(paste0(i, " : ", j))
    temp_change_col = sub_input_tree[j,2]
    if (temp_change_col == 1){
      changes_df[i,temp_change_col] = as.character(sub_input_tree[j,6])
      last_edit = as.character(sub_input_tree[j,6])
    } else {
      changes_df[i,temp_change_col] = paste0(changes_df[i,temp_change_col-1], "-", as.character(sub_input_tree[j,6]))
      last_edit = paste0(changes_df[i,temp_change_col-1], "-", as.character(sub_input_tree[j,6]))
    }
  }
  if(nrow(sub_input_tree) < ncol(changes_df)){
    for (j in (nrow(sub_input_tree) + 1): ncol(changes_df)){
      temp_change_col = j
      changes_df[i,temp_change_col] = last_edit
    }
  }
}

write.table(changes_df, paste0(output_dir,"table_of_changes_per_lineage"), sep = "\t", quote = FALSE)

#Next, let's attempt to make the tree....
attempted_cut_df = as.data.frame(matrix(data = "NONE", nrow = length(unique(input_tree$Group)), ncol = (max(input_tree$ChangeNum)) + 1), stringsAsFactors = FALSE)
rownames(attempted_cut_df) = unique(input_tree$Group)
attempted_cut_df$V1 = 1

for(i in 1:ncol(changes_df)){
  freq_col_df = as.data.frame(table(changes_df[,i]))
  freq_col_df$group_num = c(1:nrow(freq_col_df))
  
  for(j in 1:nrow(freq_col_df)){
    temp_change = as.character(freq_col_df[j,1])
    temp_group_num = as.numeric(freq_col_df[j,3])
    sub_change_df = changes_df[which(changes_df[,i] == temp_change),]
    attempted_cut_df[which(rownames(attempted_cut_df) %in% rownames(sub_change_df)),i+1] = temp_group_num
  }
}

##testing code from different script to see if we can get these numbers in order (though we may not need to??) 
## Holy shit, it just works. You write stellar, stellar code. 

just_unique_cols_file = attempted_cut_df

Modified_just_unique_cols_file = just_unique_cols_file

for (col in 2:ncol(just_unique_cols_file)){
  print(paste0("COLUMN=",col))
  last_combo = "0_1"
  max_val_in_column = 0
  for (row in 1:nrow(just_unique_cols_file)){
    #val_in_last_col = just_unique_cols_file[row,col-1]
    val_in_last_col = paste(just_unique_cols_file[row,1:(col-1)], collapse = "_")
    val_in_this_col = just_unique_cols_file[row,col]
    conc_last_this = paste0(val_in_last_col, "_", val_in_this_col)
    print(paste0(last_combo, " : ", conc_last_this, " : ", max_val_in_column))
    if (conc_last_this != last_combo){
      #print("here!")
      Modified_just_unique_cols_file[row,col] = max_val_in_column + 1
      max_val_in_column = max_val_in_column + 1
    } else {
      Modified_just_unique_cols_file[row,col] = Modified_just_unique_cols_file[row-1,col]
    }
    last_combo = conc_last_this 
  }
}

Modified_just_unique_cols_file$last_col = c(1:nrow(Modified_just_unique_cols_file))

write.table(Modified_just_unique_cols_file, paste0(output_dir, "CutDF_w_collapsed_bifurcations"), sep = "\t", quote = FALSE)





