#Modified version of 191217_CROPt_mod_MAIN_OP_for_AMB_fixes.R 
#which takes in an updated format of MAIN_OUTPUT (updated parser here)

#target = "TTAGTAGGTC"
args = commandArgs(trailingOnly=TRUE)
target = args[1]
Output_dir = args[2]
Cigar_file_name = args[3]
Main_output_file_name = args[4]

#setwd("/Users/anna/Documents/UW_GS/Shendure_Lab/Projects/Many_targets_ideas/CROPt/190722_Tree29_other7plates/TARG_plates2thru8/analysisNewPipeline/plate1")
main_op_file = read.table(paste0(Output_dir,"/",Main_output_file_name), header = FALSE, sep = "_")
#for target ordering: 
cigar_tbl = read.table(paste0(Output_dir,"/",Cigar_file_name), header = TRUE, sep = "\t")

#split V1 of main_op
temp_main_tbl = data.frame(do.call('rbind', strsplit(as.character(main_op_file$V1), ')', fixed=TRUE)))

temp_tbl = data.frame(do.call('rbind', strsplit(as.character(main_op_file$V4), '|', fixed=TRUE)))
full_tbl = as.data.frame(cbind(as.numeric(as.character(temp_main_tbl$X2)), as.character(main_op_file$V2), as.character(main_op_file$V3), as.character(temp_tbl$X1), as.character(temp_tbl$X2)))
full_tbl$cell_target = paste0(full_tbl$V2, "_", full_tbl$V3, "_", full_tbl$V4)
full_tbl$cell = paste0(full_tbl$V2, "_", full_tbl$V3)
full_tbl = full_tbl[which(full_tbl$cell %in% rownames(cigar_tbl)),]
full_tbl_targ = full_tbl[which(full_tbl$V4 == target),]
full_tbl_targ$just_edit = substr(full_tbl_targ$V5,3,length(full_tbl_targ$V5))
new_table = full_tbl_targ[,c(7,8,1)]
write.table(new_table, paste0(Output_dir,"/Table_input_to_split_target_", target), col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)

