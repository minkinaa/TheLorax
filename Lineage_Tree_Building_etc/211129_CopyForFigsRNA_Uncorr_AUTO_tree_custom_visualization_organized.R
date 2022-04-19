#This is the organized version of the script 200806_first_custom_visualizing....

library(ggplot2)
library(reshape2)
library(gridExtra)
library(cowplot)
library(dplyr)

### Potential solution to C stack issue: https://stackoverflow.com/questions/14719349/error-c-stack-usage-is-too-close-to-the-limit

#FUNCTIONS
# plot_func_whole_line = function(vect){
#   x1 = as.numeric(vect[1])
#   y1 = as.numeric(vect[2])
#   x2 = as.numeric(vect[3])
#   y2 = as.numeric(vect[4])
#   curv = curve_param*(y2-y1)
#   if (y2 > y1){
#     ang = 45
#   } else if (y2 < y1){
#     ang = 135
#   } else {
#     ang = 0
#   }
#   geom_curve(aes(x = x1, y = y1, xend = x2, yend = y2, colour = "black"), 
#              curvature = curv, ncp = 5, angle = ang, size = 2)
# }
# 
# plot_func_whole_line_SEGMENT = function(vect){
#   x1 = as.numeric(vect[1])
#   y1 = as.numeric(vect[2])
#   x2 = as.numeric(vect[3])
#   y2 = as.numeric(vect[4])
#   curv = curve_param*(y2-y1)
#   if (y2 > y1){
#     ang = 45
#   } else if (y2 < y1){
#     ang = 135
#   } else {
#     ang = 0
#   }
#   geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2, colour = "segment"))
# }

# return_fill_color = function(tbl_row, tbl_of_colors){
#   cigar_str = as.character(tbl_row[6])
#   temp_col = as.character(tbl_of_colors[which(tbl_of_colors$V1 == cigar_str),]$V2)
#   return (temp_col)
# }

#PROGRAM

args = commandArgs(trailingOnly=TRUE)
output_dir = args[1]
cut_df_file_name = args[2]
cigar_file_name = args[3]
add_extra_columns = args[4] ##should be 'false' for now
#extra_col_cigar_file = args[5] ## I think I can run w/out this argument; if not, can always run w/ cigar_file_name again
#special_labels_string = args[6] ### a string in format 1-1-1-2-3-4-5-5-5-5 etc.

###TO COMMENT OUT!! FOR TESTING!!
output_dir = "/Users/anna/Documents/UW_GS/Shendure_Lab/Projects/Many_targets_ideas/CROPt/190722_Tree29_other7plates/TARG_plates2thru8/200623_NewPipelineAnalysis/210405_ModTree42_ATAC_Analysis/211129_rna_vs_atac_figs/"
cut_df_file_name = "/Users/anna/Documents/UW_GS/Shendure_Lab/Projects/Many_targets_ideas/CROPt/190722_Tree29_other7plates/TARG_plates2thru8/200623_NewPipelineAnalysis/210326_ModifyingTree_Take2/All_Cell_Uncorrected_New43Group_cut_df_NewestGroupsAdded"
cigar_file_name = "/Users/anna/Documents/UW_GS/Shendure_Lab/Projects/Many_targets_ideas/CROPt/190722_Tree29_other7plates/TARG_plates2thru8/200623_NewPipelineAnalysis/201116_PostTreeBuild_AddAll_Cells/Plot_All_Uncorrected_Cells_W_Groups/Ordered_cigar_file.txt"
add_extra_columns = "false"
###


#output_dir = "/Users/anna/Documents/UW_GS/Shendure_Lab/Projects/Many_targets_ideas/CROPt/190722_Tree29_other7plates/TARG_plates2thru8/200623_NewPipelineAnalysis/201028_TreeAttempt1/Corrected_p1_p3/With_Added_group_splits/201110_which_new_groups_to_split/Grp6_p3H6_TCGCGTACTT-195/_M5_UnEd1_ward.D2_all_cols/"

#Step 1: Reorder cut_df & save table:
#cut_df_file_name ="/Users/anna/Documents/UW_GS/Shendure_Lab/Projects/Many_targets_ideas/CROPt/190722_Tree29_other7plates/TARG_plates2thru8/200623_NewPipelineAnalysis/201028_TreeAttempt1/Corrected_p1_p3/With_Added_group_splits/201110_which_new_groups_to_split/Grp6_p3H6_TCGCGTACTT-195/_M5_UnEd1_ward.D2_all_cols/Cells_split_into_groups_ward.D2_AtOrAbove__M5_UnEd1_ward.D2_all_cols.txt"
cut_df = read.table(cut_df_file_name, header = TRUE, sep = "\t")

alg = "ward.D2"

#PlusCounts = "PlusCounts"
PlusCounts = "NoCounts"
#suffix = paste0("AtorAbove4_M5_UnEd1_SUBTRACTED",alg,"_all_cols_", PlusCounts)
#suffix = paste0("M5_UnEd1_",alg,"_all_cols_", PlusCounts)
suffix = paste0("All_Cells_New_Groups",alg,"_all_cols_", PlusCounts)


#cigar_file = read.table("/Users/anna/Documents/UW_GS/Shendure_Lab/Projects/Many_targets_ideas/CROPt/190722_Tree29_other7plates/TARG_plates2thru8/200623_NewPipelineAnalysis/Corrected_CombinedPlates_match-5_uned-1/Combine_Like_Cells/PostDistribution_Rep_Cell_Cigar_file_AtorAbove4", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cigar_file = read.table(cigar_file_name, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

print("Dim_Cigar:")
print(dim(cigar_file))

##OPTIONAL: this adds the duplicated (and distributed) target info to the cigar file. We do this here because we don't want these targets to influence tree-building, but nice to see these for later analysis
if (add_extra_columns == "false"){
  add_extra_columns = FALSE
} else {
  add_extra_columns = TRUE
}
#add_extra_columns = FALSE
if(add_extra_columns == TRUE){
  extra_col_cigar_file = args[5]
  extra_col_cigar = read.table(extra_col_cigar_file, sep = "\t")
  names_extra_cols = colnames(extra_col_cigar)[!(colnames(extra_col_cigar) %in% colnames(cigar_file))]
  #names_extra_cols = c("TTAGTAGGTC", "TTAGTAGGTC_2", "GTGGTTGTGG", "GTGGTTGTGG_2")
  extra_col_cigar_just_new_cols = extra_col_cigar[,which(colnames(extra_col_cigar) %in% names_extra_cols)]
}
##OPTIONAL

##OPTIONAL:
make_coord_file = TRUE
if (make_coord_file == FALSE){
  df = read.table("/Users/anna/Documents/UW_GS/Shendure_Lab/Projects/Many_targets_ideas/CROPt/190722_Tree29_other7plates/TARG_plates2thru8/200623_NewPipelineAnalysis/201028_TreeAttempt1/Corrected_p1_p3/With_Added_group_splits/Tree_Coordinates_File", header = TRUE, sep = "\t")
}
##OPTIONAL

##OPTIONAL
convert_coordinates_to_square_tree_coords = TRUE
erase_points_in_cpp_script = FALSE ##TRUE = original cpp script; works well for bifurcating tree; less well for not... 

##OPTIONAL

for(i in ncol(cut_df):1){
  #print(i)
  cut_df = cut_df[order(cut_df[,i]),]
}

cut_df_file_path = paste0(output_dir, "ReorderedCutDF_", suffix)

write.table(cut_df, cut_df_file_path, sep = "\t", quote = FALSE)
write.table(cut_df, paste0(output_dir, "ReorderedCutDF_aka_cutdf"), sep = "\t", quote = FALSE)

#print("Dim CutDF:")
#print(dim(cut_df))

#Step 2: Run through cpp program (hey, why not do this from R??)
if (make_coord_file == TRUE){
  coord_file_name = paste0(output_dir, "Coordinate_file_for_tree_plot_", suffix)
  system2("touch", args = coord_file_name)
  
  if (erase_points_in_cpp_script == TRUE){
    cpp_prog_name = "/Users/anna/Documents/UW_GS/Shendure_Lab/Projects/Many_targets_ideas/CROPt/CROPt_pipeline_files/cpp_CROPt_workspace/200806_make_coordinates_for_tree_plot/Debug/200806_make_coordinates_for_tree_plot"
  } else {
    cpp_prog_name = "/Users/anna/Documents/UW_GS/Shendure_Lab/Projects/Many_targets_ideas/CROPt/CROPt_pipeline_files/cpp_CROPt_workspace/201108_make_coordinates_for_tree_plot_NO_ERASING/Debug/201108_make_coordinates_for_tree_plot_NO_ERASING"
  }
  system2(cpp_prog_name, args = c(cut_df_file_path, coord_file_name), stdout = FALSE)
  df = read.table(coord_file_name, header = TRUE, sep = "\t")
  df = df[2:nrow(df),]
}

##GENERATE TREE PLOT##
#df = read.table(coord_file_name, header = TRUE, sep = "\t")
#df = df[2:nrow(df),]
df$x1 = as.numeric(as.character(df$x1))
df$x2 = as.numeric(as.character(df$x2))
df$y1 = as.numeric(as.character(df$y1))
df$y2 = as.numeric(as.character(df$y2))


##A TEST
multiply_branches_by = 10
df$x1 = df$x1*multiply_branches_by
df$x2 = df$x2*multiply_branches_by
##A TEST

if (convert_coordinates_to_square_tree_coords == TRUE){
  df_vert = df
  df_vert[which(df_vert$y1 != df_vert$y2),]$x2 = df_vert[which(df_vert$y1 != df_vert$y2),]$x1
  df_horiz = df
  df_horiz[which(df_horiz$y1 != df_horiz$y2),]$y1 = df_horiz[which(df_horiz$y1 != df_horiz$y2),]$y2
  df = rbind(df_vert,df_horiz)
  suffix = paste0(suffix,"_SQUARE_")
}

df_ends = df[which(df$cell_name != "<NA>"),]
df_ends = df_ends[which(df_ends$y1 == df_ends$y2),]

write.table(df, paste0(output_dir, "Tree_Coordinates_File"), sep = "\t", quote = FALSE, row.names = FALSE)

#####
#Ok, here we'll attempt to add numbers to branches:
df_w_groups = df
df_w_groups$col_num = df_w_groups$x1/multiply_branches_by+2
df_w_groups$group_num = 0
df_w_groups = df_w_groups[order(-df_w_groups$y1),]
df_w_groups = df_w_groups[order(df_w_groups$x2),]
df_w_groups = df_w_groups[which(df_w_groups$x1 != df_w_groups$x2),]
df_w_groups = distinct(df_w_groups)

for(i in 2:ncol(cut_df)){
  print(i)
  temp_col = i
  vect_of_unique_groups = unique(cut_df[,i])
  #print(length(vect_of_unique_groups))
  df_w_groups[which(df_w_groups$col_num == i),]$group_num = vect_of_unique_groups
  #print(nrow(df_w_groups[which(df_w_groups$col_num == i),]))
}

df_w_groups$col_and_group = paste0(df_w_groups$col_num, ":", df_w_groups$group_num)
x_shift = 1
y_shift = .3
df_w_groups$new_x = df_w_groups$x1 + x_shift
df_w_groups$new_y = df_w_groups$y1 + y_shift

temp_vertical_df = df[which(df$x1 == df$x2),]
temp_vertical_df$x1_y1 = paste0(temp_vertical_df$x1, "_", temp_vertical_df$y1)
temp_vertical_df$x2_y2 = paste0(temp_vertical_df$x2, "_", temp_vertical_df$y2)
combined_coords = c(temp_vertical_df$x1_y1, temp_vertical_df$x2_y2)


df_w_groups_just_starts = df_w_groups
df_w_groups_just_starts$x1_y1 = paste0(df_w_groups_just_starts$x1, "_", df_w_groups_just_starts$y1)
df_w_groups_just_starts = df_w_groups_just_starts[which(df_w_groups_just_starts$x1_y1 %in% combined_coords),]
#####

####Let's also add the original groups at the ends of each branch:
df_end_points = df[which(df$x2 == max(df$x2)),]
df_end_points = distinct(df_end_points)
df_end_points = df_end_points[order(-df_end_points$y1),]
df_end_points$label = as.character(1:nrow(df_end_points))
df_end_points$x_pos = df_end_points$x1 + x_shift
df_end_points$y_pos = df_end_points$y1 + y_shift

##get list of colors
#Original colors file:
colors_file = read.table("/Users/anna/Documents/UW_GS/Shendure_Lab/Projects/Many_targets_ideas/CROPt/190722_Tree29_other7plates/TARG_plates2thru8/200623_NewPipelineAnalysis/Corrected_CombinedPlates_match-5_uned-1/COLORS_Corrected_All_Plates_AMBcorr_Xcorr_cigar_match-5_uned-1", sep = "\t", comment.char = "", stringsAsFactors = FALSE)
#Colors file with ATAC data added:
#colors_file = read.table("/Users/anna/Documents/UW_GS/Shendure_Lab/Projects/Many_targets_ideas/CROPt/191126_targATAC_full_1/201212_NewATACtargetAnalysis/AllPlates/210106_color_list_w_ATAC_added", sep = "\t", comment.char = "", stringsAsFactors = FALSE)
cols = colors_file$V2
names(cols) = colors_file$V1

##get associated cigar file
#cigar_file = read.table("/Users/anna/Documents/UW_GS/Shendure_Lab/Projects/Many_targets_ideas/CROPt/190722_Tree29_other7plates/TARG_plates2thru8/200623_NewPipelineAnalysis/Corrected_CombinedPlates_match-5_uned-1/Combine_Like_Cells/PostDistribution_Rep_Cell_Cigar_file_AtorAbove4", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
ordered_cigar_file = cigar_file[match(rownames(cut_df),rownames(cigar_file)),]
sub_ordered_cigar_file = ordered_cigar_file[which(rownames(ordered_cigar_file) %in% df_ends$cell_name),]
sub_ordered_cigar_file = sub_ordered_cigar_file[match(df_ends$cell_name, rownames(sub_ordered_cigar_file)),]

if(add_extra_columns == TRUE){
  extra_col_cigar_just_new_cols = extra_col_cigar_just_new_cols[match(rownames(ordered_cigar_file), rownames(extra_col_cigar_just_new_cols)),]
  ordered_cigar_file = cbind(ordered_cigar_file,extra_col_cigar_just_new_cols)
  ordered_cigar_file = ordered_cigar_file[,match(colnames(extra_col_cigar), colnames(ordered_cigar_file))]
}

sub_ordered_cigar_file = ordered_cigar_file[which(rownames(ordered_cigar_file) %in% df_ends$cell_name),]
sub_ordered_cigar_file = sub_ordered_cigar_file[match(df_ends$cell_name, rownames(sub_ordered_cigar_file)),]

write.table(ordered_cigar_file, paste0(output_dir, "Ordered_cigar_file.txt"), sep = "\t", quote = FALSE)

#generate new_df_ends file
df_ends$cigar_str = sub_ordered_cigar_file$TGCCTCTCGC
df_ends = df_ends[,c(1:6)]

shift_x_by = 10

whole_set_df = df_ends
last_df = df_ends
for (i in 2:ncol(sub_ordered_cigar_file)){
  print(i)
  temp_df = last_df
  temp_df$x2 = temp_df$x2 + shift_x_by
  temp_df$cigar_str = sub_ordered_cigar_file[,i]
  last_df = temp_df
  whole_set_df = rbind(whole_set_df,temp_df)
  print(nrow(whole_set_df))
}

###CHANGE THIS VALUE IF NEEDED
last_col_max = max(cut_df[,ncol(cut_df)])
last_second_to_last_col_max = max(cut_df[,(ncol(cut_df)-1)])
if (last_col_max - last_second_to_last_col_max > 1){
  cut_df_col_to_use = (ncol(cut_df) - 1)
} else {
  cut_df_col_to_use = ncol(cut_df)
}

if(length(args) > 5){
  special_labels_string = args[6]
  split_labels = as.numeric(strsplit(special_labels_string, "-")[[1]])
  df_ends = distinct(df_ends)
  df_ends = df_ends[match(rownames(cut_df), df_ends$cell_name),]
  df_ends$cut_df_values = split_labels
  cut_df$new_labels = df_ends$cut_df_values 
  cut_df_col_to_use = ncol(cut_df)
}

#cut_df_col_to_use = (ncol(cut_df) - 1)
#cut_df_col_to_use = ncol(cut_df)
cols_to_generate = max(cut_df[,cut_df_col_to_use])
#cut_df_col_to_use = 94
#cols_to_generate = cut_df_col_to_use
#cut_df_col_to_use = ncol(cut_df)

temp_cut_df_to_reorder = cut_df
temp_cut_df_to_reorder = temp_cut_df_to_reorder[match(df_ends$cell_name,rownames(temp_cut_df_to_reorder)),]
df_ends$cut_df_values = temp_cut_df_to_reorder[,cut_df_col_to_use]

write.table(df_ends, paste0(output_dir, "Coordinates_For_Tree_Endpoints_aka_df_ends"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(whole_set_df, paste0(output_dir, "Coordinates_For_Edit_Blocks_aka_whole_set_df"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(colors_file, paste0(output_dir, "Colors_file_for_edit_patterns"), sep = "\t", quote = FALSE, row.names = FALSE)

set.seed(3)
cols_groups = rainbow(cols_to_generate, s=.6, v=.9)[sample(cols_to_generate,cols_to_generate)]
names(cols_groups) = as.character(unique(df_ends$cut_df_values))


#####TEST
#df_ends$col_val = apply(df_ends, 1, return_fill_color, tbl_of_colors = colors_file)
######
curve_param = 0

dim_factor = (nrow(cut_df)/10)
#x_dim = dim_factor
y_dim = dim_factor

x_dim = ncol(cut_df)/10 + 5

if (x_dim < 7){
  x_dim = 7
}
if(y_dim < 3){
  y_dim = 3
}
print(x_dim)
print(y_dim)

add_cutdf_info = TRUE

old_whole_set_df = whole_set_df  ## to save

group_to_use =8
temp_df_ends = df_ends[which(df_ends$cut_df_values == group_to_use),]
unique(temp_df_ends$y1)
#total_number_of_cells = 50
total_number_of_cells = length(unique(temp_df_ends$y1)) 
ys_to_keep = temp_df_ends$y1[(length(temp_df_ends$y1) - total_number_of_cells + 1):length(temp_df_ends$y1)]

temp_df_ends = temp_df_ends[which(temp_df_ends$y1 %in% ys_to_keep),]

y_dim = nrow(temp_df_ends)*12

new_whole_set_df = whole_set_df[which(whole_set_df$cell_name %in% temp_df_ends$cell_name),] 
min_y = min(new_whole_set_df$y1)
new_whole_set_df$y1 = new_whole_set_df$y1 - min_y
new_whole_set_df$y2 = new_whole_set_df$y2 - min_y

test_plot = ggplot() + geom_point()
#test_plot = test_plot + apply(df, 1, plot_func_whole_line)
#test_plot = test_plot + apply(df, 1, plot_func_whole_line_SEGMENT)
#test_plot = test_plot + geom_segment(data = df, aes(x = x1, y = y1, xend = x2, yend = y2, colour = "segment"))
test_plot = test_plot + 
  theme_minimal() + theme(legend.position="none") + 
  theme(axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank()) +
  theme(panel.background = element_rect(fill = 'white', color = 'white')) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#test_plot = test_plot + geom_point(aes(x = df_ends$x2, y = df_ends$y2))
#test_plot = test_plot + geom_text(aes(x = df_ends$x2, y = df_ends$y2, label = df_ends$cell_name), hjust = 0, nudge_x = 1)
#test_plot = test_plot + geom_point(aes(x = c(df_ends$x2 + 1, df_ends$x2 + 6), y = c(df_ends$y2, df_ends$y2)), shape = 22, size = 10)
test_plot = test_plot + geom_rect(data =new_whole_set_df, mapping=aes(xmin = x2 + 1, xmax = x2+10, ymin=y2-.4, ymax=y2+.4, fill = factor(cigar_str)), color = "black", size = .05)
#test_plot = test_plot + scale_fill_manual(values = cols)

new_x_shift = shift_x_by*ncol(sub_ordered_cigar_file) + 5
#test_plot = test_plot + geom_text(aes(x = df_ends$x2 + new_x_shift + 5, y = df_ends$y2, label = df_ends$cell_name), hjust = 0, nudge_x = 8, size = 1)
#test_plot = test_plot + geom_rect(data =df_ends, mapping=aes(xmin = x2 + new_x_shift, xmax = x2 + new_x_shift+10, ymin=y2-.4, ymax=y2+.4, fill = factor(cut_df_values)), color = "grey")
#test_plot = test_plot + geom_text(aes(x = df_ends$x2 + new_x_shift, y = df_ends$y2, label = df_ends$cut_df_values), hjust = 0, nudge_x = 2, size = 1)
test_plot = test_plot + scale_fill_manual(values = c(cols,cols_groups))

#test_plot = test_plot + geom_text(aes(x = df_w_groups_just_starts$new_x, y = df_w_groups_just_starts$new_y, label = df_w_groups_just_starts$col_and_group), hjust = 0, nudge_x = 0, size = 1)
#test_plot = test_plot + geom_text(aes(x = df_end_points$x_pos, y = df_end_points$y_pos, label = df_end_points$label), hjust = 0, nudge_x = 0, size = 1)

#test_plot = test_plot + xlim(0, new_x_shift + 170)
test_plot = test_plot + xlim(0, max(new_whole_set_df$x2 + 90))

output_file_name = paste0(output_dir,"RNA_TreePlot_Group",group_to_use,"_NumCells",total_number_of_cells,".pdf")
#pdf(output_file_name, width = 43, height = 200)
#pdf(output_file_name, width = 5, height = 5)
#pdf(output_file_name, width = 5, height = 200)

pdf(output_file_name, width =x_dim/3, height = y_dim/(400))
#pdf(output_file_name, width =x_dim, height = y_dim/100)

#par(mar=c(50,50,50,50))
test_plot
dev.off()















