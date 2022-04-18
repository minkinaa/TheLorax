#install.packages("gplots")
#install.packages("gridExtra")
#install.packages("pheatmap")
library(gplots)
library(pheatmap)
library(ggplot2)
library(reshape2)
require(gridExtra)
library(gridExtra)

sample_name = "plate1"
setwd(paste0("/Users/anna/Documents/UW_GS/Shendure_Lab/Projects/Many_targets_ideas/CROPt/191126_targATAC_full_1/201212_NewATACtargetAnalysis/", sample_name))
prefix="191126_Tree29_"
#rna_file_name="190624_Tree29_plate2_combined_dedup_metrics"

rna_file_exists = FALSE
min_num_target_cutoff= 5

if (rna_file_exists == TRUE){
  rna_file_name=paste0(prefix, sample_name, "_combined_dedup_metrics_wPlate")
}

run_number = 1
clone36=TRUE
Log_file = paste0(prefix, sample_name, "_Log_", run_number, sep = "")
if(file.exists(Log_file) != FALSE)
{
  #Log_file = Log_file[-(1:length(readLines(Log_file)))]
  writeLines("", Log_file)
}
#Index_list = ""

write(paste0("Log for ", prefix, sample_name, sep = ""), Log_file)

#test_string = paste0(prefix, "_cell_by_target_correlation_outfile", sep = "")

EorU_tbl = as.data.frame(read.table(paste0(prefix, sample_name, "_cell_by_target_E_or_U_outfile", sep = ""), header = TRUE, row.names = 1, sep = "\t"))
Most_abun_umi_tbl = as.data.frame(read.table(paste0(prefix, sample_name, "_cell_by_target_most_abundant_umi_count_outfile", sep = ""), header = TRUE, row.names = 1, sep = "\t"))
Total_umi_tbl = as.data.frame(read.table(paste0(prefix, sample_name, "_cell_by_target_total_umi_count_outfile", sep = ""), header = TRUE, row.names = 1, sep = "\t"))
Cigar_tbl = as.data.frame(read.table(paste0(prefix, sample_name, "_cell_by_target_cigar_outfile", sep = ""), header = TRUE, row.names = 1, sep = "\t"))
#Corr_tbl = as.data.frame(read.table(paste0(prefix, sample_name, "_cell_by_target_correlation_outfile", sep = ""), header = TRUE, row.names = 1, sep = "\t"))

write(paste0("Total_number_of_cells_before_cutoffs=",nrow(EorU_tbl)), Log_file, append = TRUE)

if (rna_file_exists == TRUE){
  RNA_tbl = read.table(rna_file_name, header = FALSE, sep = '\t')
  RNA_tbl$V4 = log(RNA_tbl$V2, base = 2)
  hist(RNA_tbl$V4, breaks = 40, main = "Num unique reads (UMIs) per cell, sciRNAseq", xlab = "log2(Num unique UMIs per cell)")
  
  RNA_log_cutoff = 11
  write(paste0("RNA_log_cutoff=",RNA_log_cutoff), Log_file, append = TRUE)
  RNA_tbl_real = RNA_tbl[which(RNA_tbl$V4 > RNA_log_cutoff),]
  RNA_tbl_real$V5 = RNA_tbl_real$V3/RNA_tbl_real$V2
  hist(RNA_tbl_real$V5, breaks = 20)
  
  pdf(paste0(prefix, sample_name, "_RNA_duplication_rate_", run_number, ".pdf", sep = ""), width = 4, height = 4)
  hist(RNA_tbl_real$V5, breaks = 20, main = "RNA duplication rate", xlab = "Duplication rate")
  dev.off()
  
  write(paste0("Mean RNA duplication rate=",mean(RNA_tbl_real$V5)), Log_file, append = TRUE)
  
  #remove cells w/out sciRNA profiles:
  EorU_tbl_real = EorU_tbl[which(rownames(EorU_tbl) %in% RNA_tbl_real$V1),]
  Cigar_tbl_real = Cigar_tbl[which(rownames(Cigar_tbl) %in% RNA_tbl_real$V1),]
  Most_abun_umi_tbl_real = Most_abun_umi_tbl[which(rownames(Most_abun_umi_tbl) %in% RNA_tbl_real$V1),]
  Total_umi_tbl_real = Total_umi_tbl[which(rownames(Total_umi_tbl) %in% RNA_tbl_real$V1),]
  
  write(paste0("Total_number_of_cells_after_RNA_cutoff=",nrow(EorU_tbl_real)), Log_file, append = TRUE)
  
} else {
  EorU_tbl_real = EorU_tbl
  Cigar_tbl_real = Cigar_tbl
  Most_abun_umi_tbl_real = Most_abun_umi_tbl
  Total_umi_tbl_real = Total_umi_tbl
  
  write(paste0("Total_number_of_cells_after_RNA_cutoff=",nrow(EorU_tbl_real)), Log_file, append = TRUE)

}

#remove any cells w/ low number of targets captured:
EorU_tbl_real$num_non_Xs = 0

for (row in 1:nrow(EorU_tbl_real))
{
  EorU_tbl_real[row,ncol(EorU_tbl_real)] = ncol(EorU_tbl_real) - sum(EorU_tbl_real[row,] == "X")
}

hist(EorU_tbl_real$num_non_Xs, breaks = 20, main = "Number of targets captured per cell")

pdf(paste0(prefix, sample_name, "_histogram_num_targets_per_cell_", run_number, ".pdf", sep = ""), width = 6, height = 4)
hist(EorU_tbl_real$num_non_Xs, breaks = 20, main = "Number of targets captured per cell", xlab = "Number of targets out of 45")
dev.off()

##test
#EorU_tbl_real$SUM = EorU_tbl_real$num_A + EorU_tbl_real$num_E + EorU_tbl_real$num_U
#Total_umi_tbl_real$SUM = rowSums(Total_umi_tbl_real)
#Total_umi_tbl_real$mean_UMIs_per_target = Total_umi_tbl_real$SUM/EorU_tbl_real$num_non_Xs
#plot(Total_umi_tbl_real$mean_UMIs_per_target, EorU_tbl_real$num_non_Xs, xlab = "Mean_UMIs_per_target", ylab = "Number of targets")
#Total_umi_tbl_real$real = FALSE
#rna_real_actual = RNA_tbl[which(RNA_tbl$V4 > RNA_log_cutoff),]
#Total_umi_tbl_real[which(rownames(Total_umi_tbl_real) %in% rna_real_actual$V1),]$real = TRUE
#ggplot(Total_umi_tbl_real, aes(mean_UMIs_per_target, num_targets_captured)) + geom_point(aes(colour = factor(real)))
#RNA_tbl$has_rna = "UNKNOWN"
#RNA_tbl[which(RNA_tbl$V1 %in% rownames(Total_umi_tbl_real[which(Total_umi_tbl_real$real == TRUE),])),]$has_rna = "YES"
#RNA_tbl[which(RNA_tbl$V1 %in% rownames(Total_umi_tbl_real[which(Total_umi_tbl_real$real == FALSE),])),]$has_rna = "NO"
#RNA_tbl_Y_N = RNA_tbl[(which(RNA_tbl$has_rna != "UNKNOWN")),]
#ggplot(RNA_tbl_Y_N[which(RNA_tbl_Y_N$V1 %in% rownames(Total_umi_tbl_real[which(Total_umi_tbl_real$num_targets_captured > 19),])),], aes(has_rna,V2)) + geom_boxplot()
#Total_umi_tbl_real$RT_index = substr(rownames(Total_umi_tbl_real), nchar(rownames(Total_umi_tbl_real))- 9, nchar(rownames(Total_umi_tbl_real)))
#temp_tbl = Total_umi_tbl_real[which(Total_umi_tbl_real$num_targets_captured > 19 & Total_umi_tbl_real$real == FALSE),]
#table(temp_tbl$RT_index)
##test
#min_num_target_cutoff= 10
write(paste0("min_num_target_cutoff=",min_num_target_cutoff), Log_file, append = TRUE)
EorU_tbl_real = EorU_tbl_real[which(EorU_tbl_real$num_non_Xs > min_num_target_cutoff),]
#EorU_tbl_real= EorU_tbl_real[,-which(colnames(EorU_tbl_real) %in% "CACTTCGAGG")]
if (clone36 == TRUE)
{
  EorU_tbl_real= EorU_tbl_real[,-which(colnames(EorU_tbl_real) %in% "CACTTCGAGG")]
}
Cigar_tbl_real = Cigar_tbl_real[which(rownames(Cigar_tbl_real) %in% rownames(EorU_tbl_real)), colnames(Cigar_tbl_real) %in% colnames(EorU_tbl_real)]
Most_abun_umi_tbl_real = Most_abun_umi_tbl_real[which(rownames(Most_abun_umi_tbl_real) %in% rownames(EorU_tbl_real)), colnames(Most_abun_umi_tbl_real) %in% colnames(EorU_tbl_real)]
Total_umi_tbl_real = Total_umi_tbl_real[which(rownames(Total_umi_tbl_real) %in% rownames(EorU_tbl_real)),colnames(Total_umi_tbl_real) %in% colnames(EorU_tbl_real)]

EorU_tbl_real = EorU_tbl_real[,1:ncol(EorU_tbl_real) - 1]

write(paste0("Total_number_of_cells_after_min_target_cutoff=",nrow(EorU_tbl_real)), Log_file, append = TRUE)

#Fix Cigar tbl to include AMB for ambiguous plots:
Cigar_tbl_real[]= lapply(Cigar_tbl_real, as.character)

for (i in 1:nrow(Cigar_tbl_real))
{
  for (j in 1:ncol(Cigar_tbl_real))
  {
    if (EorU_tbl_real[i,j] == "A")
    {
      Cigar_tbl_real[i,j] = "AMB"
    }
  }
}

#reorder targets by fraction of cells in which they are expressed:
EorU_tbl_real_temp = EorU_tbl_real
EorU_tbl_real = EorU_tbl_real[,order(-colSums(EorU_tbl_real_temp != "X"))]
Cigar_tbl_real = Cigar_tbl_real[,order(-colSums(EorU_tbl_real_temp != "X"))]
Most_abun_umi_tbl_real = Most_abun_umi_tbl_real[,order(-colSums(EorU_tbl_real_temp != "X"))]
Total_umi_tbl_real=Total_umi_tbl_real[,order(-colSums(EorU_tbl_real_temp != "X"))]

#add another matrix that's EorU, but includes insertion/deletion information instead of the E's:
Del_or_Ins_EorU_tbl_real = EorU_tbl_real
Del_or_Ins_EorU_tbl_real[] = lapply(Del_or_Ins_EorU_tbl_real, as.character)

for (i in 1:nrow(Cigar_tbl_real))
{
  for (j in 1:ncol(Cigar_tbl_real))
  {
    if (grepl("D", Cigar_tbl_real[i,j]) == TRUE)
    {
      Del_or_Ins_EorU_tbl_real[i,j] = "D"
    }
    else if (grepl("I", Cigar_tbl_real[i,j]) == TRUE)
    {
      Del_or_Ins_EorU_tbl_real[i,j] = "I"
    }
    else if (grepl("N", Cigar_tbl_real[i,j]) == TRUE)
    {
      Del_or_Ins_EorU_tbl_real[i,j] = "N"
    }
  }
}
#Del_or_Ins_EorU_tbl_real[] = lapply(Del_or_Ins_EorU_tbl_real, as.factor)

#plot insertion lengths
Ins_length_vect = vector()

for(i in 1:nrow(Del_or_Ins_EorU_tbl_real))
{
  for(j in 1:ncol(Del_or_Ins_EorU_tbl_real))
  {
    if(Del_or_Ins_EorU_tbl_real[i,j] == "I")
    {
      length_of_insertion = nchar(strsplit(Cigar_tbl_real[i,j],":")[[1]][2])
      if (length_of_insertion > 60)
      {
        print(Cigar_tbl_real[i,j])
      }
      Ins_length_vect = c(Ins_length_vect,length_of_insertion)
    }
  }
}

pdf(paste0(prefix, sample_name, "_histogram_of_insertion_lengths_", run_number, ".pdf", sep = ""), width = 6, height = 4)
hist(Ins_length_vect, breaks = 50, main = "Insertion length distribution", xlab = "Length of Insertion")
dev.off()

Del_or_Ins_EorU_tbl_real[] = lapply(Del_or_Ins_EorU_tbl_real, as.factor)

#make cell by EorU plot:
EorU_as_matrix = as.matrix(EorU_tbl_real)
EorU_as_matrix_melt = melt(EorU_as_matrix)
EorU_as_matrix_melt$value = as.factor(EorU_as_matrix_melt$value)
EorU_as_matrix_melt$width = 1
EorU_as_matrix_melt$height = .8

plot_EorU_tbl = ggplot(EorU_as_matrix_melt, aes(x = Var2, y = Var1)) +
  geom_tile(aes(fill=value, width = EorU_as_matrix_melt$width, height = EorU_as_matrix_melt$height))
cols = c("X" = "black", "E" = "RED", "U" = "yellow", "A" = "blue")

plot_EorU_tbl = plot_EorU_tbl + scale_fill_manual(values = cols) + theme_minimal() + theme(legend.position="right") + theme(axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x=element_text(angle=90,hjust=1), axis.text.y=element_text(size = 6))

pdf(paste0(prefix, sample_name, "_EorU_matrix_plot_", run_number, ".pdf", sep = ""), width = 10, height = nrow(EorU_tbl_real)/70*5)
plot_EorU_tbl
dev.off()

pdf(paste0(prefix, sample_name, "_EorU_matrix_plot_condensed_", run_number, ".pdf", sep = ""), width = 10, height = 10)
plot_EorU_tbl
dev.off()

#count total number of cells in which each target is captured:
pdf(paste0(prefix, sample_name, "_EorU_captured_target_bargraphs_", run_number, ".pdf", sep = ""))
barplot(colSums(EorU_tbl_real != "X")/colSums(EorU_tbl_real != "something"), las=2, xaxt='n', main = "Fraction of cells in which target captured")
barplot((colSums(EorU_tbl_real == "E") + colSums(EorU_tbl_real == "A"))/colSums(EorU_tbl_real != "X"), las=2, xaxt='n', main = "Fraction of captured targets which are edited or ambiguous")
barplot((colSums(EorU_tbl_real == "E"))/colSums(EorU_tbl_real != "X"), las=2, xaxt='n', main = "Fraction of captured targets which are edited")
barplot(colSums(EorU_tbl_real == "A")/colSums(EorU_tbl_real != "X"), las=2, main = "Fraction of captured targets which are ambiguous")
dev.off()

#make cell by total UMIs plot:
Total_umi_tbl_real_wNAs = Total_umi_tbl_real
Total_umi_tbl_real_wNAs[Total_umi_tbl_real_wNAs == 0] = NA
Total_UMIs_as_matrix = as.matrix(Total_umi_tbl_real_wNAs)
Total_UMIs_as_matrix_melt = melt(Total_UMIs_as_matrix)
Total_UMIs_as_matrix_melt$width = 1
Total_UMIs_as_matrix_melt$height = .8

plot_Total_UMIs_tbl = ggplot(Total_UMIs_as_matrix_melt, aes(x = Var2, y = Var1)) +
  geom_tile(aes(fill=value, width = Total_UMIs_as_matrix_melt$width, height = Total_UMIs_as_matrix_melt$height)) + scale_fill_gradientn(colours = c("red","yellow","white"), na.value = "black")

plot_Total_UMIs_tbl = plot_Total_UMIs_tbl + theme_minimal() + theme(legend.position="right") + theme(axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x=element_text(angle=90,hjust=1), axis.text.y=element_text(size = 6))

pdf(paste0(prefix, sample_name, "_Total_UMIs_matrix_plot_", run_number, ".pdf", sep = ""), width = 10, height = nrow(EorU_tbl_real)/70*5)
plot_Total_UMIs_tbl
dev.off()

pdf(paste0(prefix, sample_name, "_Total_UMIs_matrix_plot_condensed_", run_number, ".pdf", sep = ""), width = 10, height = 10)
plot_Total_UMIs_tbl
dev.off()

#barplot(colSums(Total_umi_tbl_real)/colSums(EorU_tbl_real != "X"))
#count total number of UMIs per target:
pdf(paste0(prefix, sample_name, "_total_UMIs_captured_per_target_bargraphs_", run_number, ".pdf", sep = ""), width = 5, height = 10)
par(oma=c(8,.5,2,.5))
test = layout(matrix(c(1,2,3,4,5,6),ncol=1), widths=c(4,4,4,4,4,4), heights=c(2,1,1,1,1,1,2), TRUE)
par(mai = c(.3, .5, .3, .5))
boxplot(Total_umi_tbl_real, las=2, xaxt='n', main = "UMIs captured per cell (per target)")
par(mai = c(.3, .5, .3, .5))
barplot(colSums(Total_umi_tbl_real), las=2, xaxt='n', main = "Total UMIs captured")
barplot(colSums(Total_umi_tbl_real)/colSums(EorU_tbl_real != "X"), las=2, xaxt='n', main = "Mean total UMIs per captured target")
Total_umi_tbl_real_temp = Total_umi_tbl_real
Total_umi_tbl_real_temp[EorU_tbl_real != "E"] = 0
barplot(colSums(Total_umi_tbl_real_temp)/colSums(EorU_tbl_real == "E"), las=2, xaxt='n', main = "Mean total UMIs per edited target")
Total_umi_tbl_real_temp = Total_umi_tbl_real
Total_umi_tbl_real_temp[EorU_tbl_real != "A"] = 0
barplot(colSums(Total_umi_tbl_real_temp)/colSums(EorU_tbl_real == "A"), las=2, xaxt='n', main = "Mean total UMIs per ambiguous target")
Total_umi_tbl_real_temp = Total_umi_tbl_real
Total_umi_tbl_real_temp[EorU_tbl_real != "U"] = 0
barplot(colSums(Total_umi_tbl_real_temp)/colSums(EorU_tbl_real == "U"), las=2, main = "Mean total UMIs per unedited target")
dev.off()

#make cell by most abundant UMIs plot:
Most_abun_umi_tbl_real_wNAs = Most_abun_umi_tbl_real
Most_abun_umi_tbl_real_wNAs[Most_abun_umi_tbl_real_wNAs == 0] = NA
Most_abun_UMIs_as_matrix = as.matrix(Most_abun_umi_tbl_real_wNAs)
Most_abun_UMIs_as_matrix_melt = melt(Most_abun_UMIs_as_matrix)
Most_abun_UMIs_as_matrix_melt$width = 1
Most_abun_UMIs_as_matrix_melt$height = .8

plot_Most_abun_UMIs_tbl = ggplot(Most_abun_UMIs_as_matrix_melt, aes(x = Var2, y = Var1)) +
  geom_tile(aes(fill=value, width = Most_abun_UMIs_as_matrix_melt$width, height = Most_abun_UMIs_as_matrix_melt$height)) + scale_fill_gradientn(colours = c("red","yellow","white"), na.value = "black")

plot_Most_abun_UMIs_tbl = plot_Most_abun_UMIs_tbl + theme_minimal() + theme(legend.position="right") + theme(axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x=element_text(angle=90,hjust=1), axis.text.y=element_text(size = 6))

pdf(paste0(prefix, sample_name, "_Most_abun_UMIs_matrix_plot_", run_number, ".pdf", sep = ""), width = 10, height = nrow(EorU_tbl_real)/70*5)
plot_Most_abun_UMIs_tbl
dev.off()

pdf(paste0(prefix, sample_name, "_Most_abun_UMIs_matrix_plot_condensed_", run_number, ".pdf", sep = ""), width = 10, height = 10)
plot_Most_abun_UMIs_tbl
dev.off()

#make cell by most abundant UMIs plot w/ log scale
Most_abun_UMIs_as_matrix_melt$log2 = log(Most_abun_UMIs_as_matrix_melt$value, base = 2)
plot_Most_abun_UMIs_tbl = ggplot(Most_abun_UMIs_as_matrix_melt, aes(x = Var2, y = Var1)) +
  geom_tile(aes(fill=log2, width = Most_abun_UMIs_as_matrix_melt$width, height = Most_abun_UMIs_as_matrix_melt$height)) + scale_fill_gradientn(colours = c("red","yellow","white"), na.value = "black")

plot_Most_abun_UMIs_tbl = plot_Most_abun_UMIs_tbl + theme_minimal() + theme(legend.position="right") + theme(axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x=element_text(angle=90,hjust=1), axis.text.y=element_text(size = 6))

pdf(paste0(prefix, sample_name, "_Most_abun_UMIs_matrix_plot_LOG2_", run_number, ".pdf", sep = ""), width = 10, height = nrow(EorU_tbl_real)/70*5)
plot_Most_abun_UMIs_tbl
dev.off()

pdf(paste0(prefix, sample_name, "_Most_abun_UMIs_matrix_plot_condensed_LOG2_", run_number, ".pdf", sep = ""), width = 10, height = 10)
plot_Most_abun_UMIs_tbl
dev.off()

#plot fractions of edited/unedited/etc.targets as color barplot
EorU_tbl_real$num_E = 0
EorU_tbl_real$num_A = 0
EorU_tbl_real$num_U = 0
EorU_tbl_real$num_X = 0

for (row in 1:nrow(EorU_tbl_real))
{
  EorU_tbl_real[row,which(colnames(EorU_tbl_real) %in% "num_X")] = sum(EorU_tbl_real[row,] == "X")
  EorU_tbl_real[row,which(colnames(EorU_tbl_real) %in% "num_E")] = sum(EorU_tbl_real[row,] == "E")
  EorU_tbl_real[row,which(colnames(EorU_tbl_real) %in% "num_A")] = sum(EorU_tbl_real[row,] == "A")
  EorU_tbl_real[row,which(colnames(EorU_tbl_real) %in% "num_U")] = sum(EorU_tbl_real[row,] == "U")
}

EorU_tbl_real$EA = EorU_tbl_real$num_A + EorU_tbl_real$num_E
EorU_tbl_real_ordered_by_edits = EorU_tbl_real[order(-EorU_tbl_real$EA),]
fraction_unedited_cells = nrow(EorU_tbl_real_ordered_by_edits[which(EorU_tbl_real_ordered_by_edits$EA == 0),])/nrow(EorU_tbl_real_ordered_by_edits)

write(paste0("Fraction of cells w/ zero edits =",fraction_unedited_cells), Log_file, append = TRUE)

EorU_tbl_real_ordered_by_edits = as.matrix(EorU_tbl_real_ordered_by_edits)
EorU_tbl_real_ordered_by_edits_melt = melt(EorU_tbl_real_ordered_by_edits)
EorU_tbl_real_ordered_by_edits_melt$value = as.factor(EorU_tbl_real_ordered_by_edits_melt$value)
EorU_tbl_real_ordered_by_edits_melt$width = 1
EorU_tbl_real_ordered_by_edits_melt$height = .8

EorU_ordered_melt = rbind(EorU_tbl_real_ordered_by_edits_melt[which(EorU_tbl_real_ordered_by_edits_melt$value == "X"),], EorU_tbl_real_ordered_by_edits_melt[which(EorU_tbl_real_ordered_by_edits_melt$value == "A"),], EorU_tbl_real_ordered_by_edits_melt[which(EorU_tbl_real_ordered_by_edits_melt$value == "U"),], EorU_tbl_real_ordered_by_edits_melt[which(EorU_tbl_real_ordered_by_edits_melt$value == "E"),])

EorU_ordered_melt$value = factor(EorU_ordered_melt$value, levels = c("X", "U", "A", "E"))
EorU_ordered_barplot = ggplot(EorU_ordered_melt, aes(Var1)) + geom_bar(aes(fill=value))
cols = c("X" = "black", "E" = "RED", "U" = "yellow", "A" = "blue")
EorU_ordered_barplot = EorU_ordered_barplot + scale_fill_manual(values = cols) + theme_minimal() + theme(legend.position="none") + theme(axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x=element_text(angle=90,hjust=1, size = 5))
pdf(paste0(prefix, sample_name, "_AUEX_barplot_", run_number, ".pdf", sep = ""), width = nrow(EorU_tbl_real)/50*5, height = 10)
EorU_ordered_barplot
dev.off()

#make EorU matrix plot w/ rows ordered by edits:
EorU_tbl_real_ordered = as.data.frame(EorU_tbl_real_ordered_by_edits[,c(1:(ncol(EorU_tbl_real_ordered_by_edits) - 5))])
EorU_tbl_real_ordered_mat = as.matrix(EorU_tbl_real_ordered)
EorU_tbl_real_ordered_mat_melt = melt(EorU_tbl_real_ordered_mat)
EorU_tbl_real_ordered_mat_melt$value = as.factor(EorU_tbl_real_ordered_mat_melt$value)
EorU_tbl_real_ordered_mat_melt$width = 1
EorU_tbl_real_ordered_mat_melt$height = .8

plot_EorU_tbl_ordered = ggplot(EorU_tbl_real_ordered_mat_melt, aes(x = Var2, y = Var1)) +
  geom_tile(aes(fill=value, width = EorU_tbl_real_ordered_mat_melt$width, height = EorU_tbl_real_ordered_mat_melt$height))
cols = c("X" = "black", "E" = "RED", "U" = "yellow", "A" = "blue")

plot_EorU_tbl_ordered = plot_EorU_tbl_ordered + scale_fill_manual(values = cols) + theme_minimal() + theme(legend.position="right") + theme(axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x=element_text(angle=90,hjust=1), axis.text.y=element_text(size = 6))

pdf(paste0(prefix, sample_name, "_EorU_matrix_OREDERD_by_EA_plot_", run_number, ".pdf", sep = ""), width = 10, height = nrow(EorU_tbl_real)/70*5)
plot_EorU_tbl_ordered
dev.off()

#plot fractions of del/ins/unedited/etc.targets as color barplot
Del_or_Ins_EorU_tbl_real$num_D = 0
Del_or_Ins_EorU_tbl_real$num_I = 0
Del_or_Ins_EorU_tbl_real$num_N = 0
Del_or_Ins_EorU_tbl_real$num_A = 0
Del_or_Ins_EorU_tbl_real$num_U = 0
Del_or_Ins_EorU_tbl_real$num_X = 0

for (row in 1:nrow(Del_or_Ins_EorU_tbl_real))
{
  Del_or_Ins_EorU_tbl_real[row, which(colnames(Del_or_Ins_EorU_tbl_real) %in% "num_D")] = sum(Del_or_Ins_EorU_tbl_real[row,] == "D")
  Del_or_Ins_EorU_tbl_real[row, which(colnames(Del_or_Ins_EorU_tbl_real) %in% "num_I")] = sum(Del_or_Ins_EorU_tbl_real[row,] == "I")
  Del_or_Ins_EorU_tbl_real[row, which(colnames(Del_or_Ins_EorU_tbl_real) %in% "num_N")] = sum(Del_or_Ins_EorU_tbl_real[row,] == "N")
  Del_or_Ins_EorU_tbl_real[row, which(colnames(Del_or_Ins_EorU_tbl_real) %in% "num_A")] = sum(Del_or_Ins_EorU_tbl_real[row,] == "A")
  Del_or_Ins_EorU_tbl_real[row, which(colnames(Del_or_Ins_EorU_tbl_real) %in% "num_U")] = sum(Del_or_Ins_EorU_tbl_real[row,] == "U")
  Del_or_Ins_EorU_tbl_real[row, which(colnames(Del_or_Ins_EorU_tbl_real) %in% "num_X")] = sum(Del_or_Ins_EorU_tbl_real[row,] == "X")
}

Del_or_Ins_EorU_tbl_real$DINA = Del_or_Ins_EorU_tbl_real$num_D + Del_or_Ins_EorU_tbl_real$num_I + Del_or_Ins_EorU_tbl_real$num_N + Del_or_Ins_EorU_tbl_real$num_A
Del_or_Ins_EorU_tbl_real_ordered = Del_or_Ins_EorU_tbl_real[order(-Del_or_Ins_EorU_tbl_real$DINA),]

Del_or_Ins_EorU_tbl_real_ordered = as.matrix(Del_or_Ins_EorU_tbl_real_ordered)
Del_or_Ins_EorU_tbl_real_ordered_melt = melt(Del_or_Ins_EorU_tbl_real_ordered)
Del_or_Ins_EorU_tbl_real_ordered_melt$value = as.factor(Del_or_Ins_EorU_tbl_real_ordered_melt$value)
Del_or_Ins_EorU_tbl_real_ordered_melt$width = 1
Del_or_Ins_EorU_tbl_real_ordered_melt$height = .8

Del_or_Ins_EorU_ordered_melt = Del_or_Ins_EorU_tbl_real_ordered_melt[which(Del_or_Ins_EorU_tbl_real_ordered_melt$value %in% c("A", "D","I", "N", "U", "X")),]
#Del_or_Ins_EorU_ordered_melt = rbind() ##may need to add this??
Del_or_Ins_EorU_ordered_melt$value = factor(Del_or_Ins_EorU_ordered_melt$value, levels = c("X","U", "A","N","I","D"))
Del_or_Ins_EorU_ordered_barplot = ggplot(Del_or_Ins_EorU_ordered_melt, aes(Var1)) + geom_bar(aes(fill=value))
cols = c("X" = "black", "D" = "RED", "U" = "yellow", "A" = "blue", "I" = "violet", "N" = "orange")
Del_or_Ins_EorU_ordered_barplot = Del_or_Ins_EorU_ordered_barplot + scale_fill_manual(values = cols) + theme_minimal() + theme(legend.position="none") + theme(axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x=element_text(angle=90,hjust=1, size = 5))
pdf(paste0(prefix, sample_name, "_AUDINX_barplot_", run_number, ".pdf", sep = ""), width = nrow(Del_or_Ins_EorU_tbl_real)/50*5, height = 10)
Del_or_Ins_EorU_ordered_barplot
dev.off()

#make AUDINX matrix plot w/ rows ordered by edits:
Del_or_Ins_EorU_tbl_real_ordered = as.data.frame(Del_or_Ins_EorU_tbl_real_ordered[,c(1:(ncol(Del_or_Ins_EorU_tbl_real_ordered) - 7))])
Del_or_Ins_EorU_tbl_real_ordered_mat = as.matrix(Del_or_Ins_EorU_tbl_real_ordered)
Del_or_Ins_EorU_tbl_real_ordered_mat_melt = melt(Del_or_Ins_EorU_tbl_real_ordered_mat)
Del_or_Ins_EorU_tbl_real_ordered_mat_melt$value = as.factor(Del_or_Ins_EorU_tbl_real_ordered_mat_melt$value)
Del_or_Ins_EorU_tbl_real_ordered_mat_melt$width = 1
Del_or_Ins_EorU_tbl_real_ordered_mat_melt$height = .8

plot_Del_or_Ins_EorU_ordered = ggplot(Del_or_Ins_EorU_tbl_real_ordered_mat_melt, aes(x = Var2, y = Var1)) +
  geom_tile(aes(fill=value, width = Del_or_Ins_EorU_tbl_real_ordered_mat_melt$width, height = Del_or_Ins_EorU_tbl_real_ordered_mat_melt$height))
cols = c("X" = "black", "D" = "RED", "U" = "yellow", "A" = "blue", "I" = "violet", "N" = "orange")

plot_Del_or_Ins_EorU_ordered = plot_Del_or_Ins_EorU_ordered + scale_fill_manual(values = cols) + theme_minimal() + theme(legend.position="right") + theme(axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x=element_text(angle=90,hjust=1), axis.text.y=element_text(size = 6))

pdf(paste0(prefix, sample_name, "_AUDINX_matrix_OREDERD_by_ADIN_plot_", run_number, ".pdf", sep = ""), width = 10, height = nrow(Del_or_Ins_EorU_tbl_real)/70*5)
plot_Del_or_Ins_EorU_ordered
dev.off()

write.table(EorU_tbl_real, file = paste0("Table_", prefix, sample_name, "_EorU_tbl_real_", run_number, sep = ""), quote = FALSE, sep = "\t")
write.table(Cigar_tbl_real, file = paste0("Table_", prefix, sample_name, "_Cigar_tbl_real_", run_number, sep = ""), quote = FALSE, sep = "\t")
write.table(Del_or_Ins_EorU_tbl_real_ordered, file = paste0("Table_", prefix, sample_name, "_Del_or_Ins_EorU_tbl_real_ordered_", run_number, sep = ""), quote = FALSE, sep = "\t")

EorU_tbl_real$num_edited = 0
for (i in 1:nrow(EorU_tbl_real))
{
  EorU_tbl_real[i,ncol(EorU_tbl_real)] = sum(EorU_tbl_real[i,] == "E") + sum(EorU_tbl_real[i,] == "A")
}

hist(EorU_tbl_real$num_edited, breaks = 20)

#subset singlets
Singlet_EorU_tbl_real = EorU_tbl_real[which(EorU_tbl_real$num_A < EorU_tbl_real$num_E),]

write.table(Singlet_EorU_tbl_real, file = paste0("Table_", prefix, sample_name, "_Singlet_EorU_tbl_real_", run_number, sep = ""), quote = FALSE, sep = "\t")

###Just continue; here we subset to cells which have at least one edit
min_number_of_edits = 1
greater_than_edits = min_number_of_edits-1

Edited_singlet = Singlet_EorU_tbl_real[which(Singlet_EorU_tbl_real$num_edited > greater_than_edits),]
Edited_Cigar_tbl_real = Cigar_tbl_real[which(rownames(Cigar_tbl_real) %in% rownames(Edited_singlet)),]

library(plyr)
col_sums = as.data.frame(ldply(Edited_Cigar_tbl_real, function(c) sum(c=="X")))
col_sums$fraction_X = col_sums$V1/nrow(Edited_Cigar_tbl_real)
hist(col_sums$fraction_X, breaks = 50)
barplot(col_sums$fraction_X)
##check to make sure non between .5 & 1 -- if some are, change the .5 to something higher)
col_sums$fraction_X

only_less_than_50 = col_sums[which(col_sums$fraction_X < .8),]

only_colnames = colnames(Edited_Cigar_tbl_real)[colnames(Edited_Cigar_tbl_real) %in% only_less_than_50$.id]

Exclude_high_X_tbl = Edited_Cigar_tbl_real[,which(colnames(Edited_Cigar_tbl_real) %in% only_colnames)]

Edited_Cigar_tbl_real = Exclude_high_X_tbl
temp_cigar_tbl_readl = Cigar_tbl_real
#Cigar_tbl_real =temp_cigar_tbl_readl
Cigar_tbl_real = Edited_Cigar_tbl_real

write.table(Cigar_tbl_real, file = paste0("Cigar_tbl_for_corr_input_greater_than_", greater_than_edits,"_edited"), col.names = FALSE, quote = FALSE, sep = "\t")
write.table(Cigar_tbl_real, file = paste0("Cigar_tbl_for_corr_input_greater_than_", greater_than_edits,"_edited_wCOLNAMES"), col.names = TRUE, quote = FALSE, sep = "\t")
write.table(cbind(colnames(Cigar_tbl_real)), file = paste0("Cigar_tbl_for_corr_input_greater_than_", greater_than_edits,"_edited_COLNAMES"), col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
write.table(rownames(Cigar_tbl_real), file = paste0("Cigar_tbl_for_corr_input_greater_than_", greater_than_edits,"_edited_ROWNAMES"), col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

##reorder columns based on past stuff;
tbl_w_cols_in_order = read.table("/Users/anna/Documents/UW_GS/Shendure_Lab/Projects/Many_targets_ideas/CROPt/191120_first_targ_atac/analysis/CIGAR_ordered_columns_w_names", header = TRUE, sep = "\t")
Cigar_tbl_real_ordered_cols = Cigar_tbl_real[,match(colnames(tbl_w_cols_in_order), colnames(Cigar_tbl_real))]
write.table(Cigar_tbl_real_ordered_cols, file = paste0("Cigar_tbl_for_corr_input_greater_than_", greater_than_edits,"_edited_ORDEREDwCOLNAMES"), col.names = TRUE, quote = FALSE, sep = "\t")


##### 06/24/2020 Note: I think this is as far as we want to go! The rest is grouping/tree-building and is done later... 









#### HERE RUN CPP SCRIPT 
# ** cd into output folder!
# check that parameters are what you want them to be!
# /Users/anna/Documents/UW_GS/Shendure_Lab/Projects/Many_targets_ideas/CROPt/CROPt_pipeline_files/Scripts/1900820_distance_matrix_bash_wrapper_script_2.txt Cigar_tbl_for_corr_input_greater_than_0_edited

Run_num = 1 ##Have this match cpp script

normalized_dist_mat = read.table("Non_Normalized_Distance_Mat_1", header = FALSE, sep = "\t")
rownames(normalized_dist_mat) = normalized_dist_mat$V1
normalized_dist_mat = normalized_dist_mat[,c(2:ncol(normalized_dist_mat))]
colnames(normalized_dist_mat) = rownames(normalized_dist_mat)

#plot = pheatmap(normalized_dist_mat, fontsize = 1, silent = TRUE)
plot = pheatmap(normalized_dist_mat, clustering_method = "ward.D", fontsize = 1, silent = TRUE)
row_names_in_order = as.data.frame(rownames(normalized_dist_mat[plot$tree_row[["order"]],]))
colnames(row_names_in_order) = c("cells_in_order")
rownames(row_names_in_order) = row_names_in_order$cells_in_order

#pdf(paste0(Output_dir,"Heatmap_NON_normalized_Run_num", Run_num), width = 10, height = 10)
pdf(paste0("Heatmap_NON_normalized_Run_num", Run_num, "_wardD"), width = 10, height = 10)
plot
dev.off()

#make the tree
cut_tree_1 = cutree(plot$tree_row, k = 1)
cut_tree_df = as.data.frame(cut_tree_1)

max_num_groups = nrow(normalized_dist_mat)

for(i in 2:max_num_groups)
{
  #col_name = paste0("cut_tree_", as.character(i))
  temp =  cutree(plot$tree_row, k = i)
  cut_tree_df = cbind(cut_tree_df, temp)
}
colnames(cut_tree_df) = c(1:max_num_groups)

write.table(cut_tree_df, file = paste0("cells_split_into_groups_wardD_Run_Num_", Run_num), quote = FALSE, sep = "\t")

###FROM HERE, go to script: 191102_make_LG_group_plots_JUST....




###########STOP###############














######

#test:



######

#



####STOP HERE FOR NOW~~~~~

##Add a step here to read the tables saved above & change values to characer, numertic, whatever is needed.

#generate file of edited cells for corr matrix: ***DO THIS AT THE END!!!!!

### If not doing immediately, start here:
setwd("/Users/anna/Documents/UW_GS/Shendure_Lab/Projects/Many_targets_ideas/CROPt/191120_first_targ_atac/analysis")
sample_name = "Tree29all"
prefix="191120_Tree29_"
run_number = 1
EorU_tbl_real = as.data.frame(read.table(paste0("Table_", prefix, sample_name, "_EorU_tbl_real_", run_number, sep = "")))
Cigar_tbl_real = as.data.frame(read.table(paste0("Table_", prefix, sample_name, "_Cigar_tbl_real_", run_number, sep = "")))
Del_or_Ins_EorU_tbl_real_ordered = as.data.frame(paste0("Table_", prefix, sample_name, "_Del_or_Ins_EorU_tbl_real_ordered_", run_number, sep = ""))



####Or, could just upload Singlet table
Singlet_EorU_tbl_real = as.data.frame(read.table("/Users/anna/Documents/UW_GS/Shendure_Lab/Projects/Many_targets_ideas/CROPt/191126_targATAC_full_1/analysis/plate1/Table_191126_Tree29_plate1_Singlet_EorU_tbl_real_1"))

###

#EorU_tbl_real$num_edited = 0
#for (i in 1:nrow(EorU_tbl_real))
#{
#  EorU_tbl_real[i,ncol(EorU_tbl_real)] = sum(EorU_tbl_real[i,] == "E") + sum(EorU_tbl_real[i,] == "A")
#}

#hist(EorU_tbl_real$num_edited, breaks = 20)
#Singlet_EorU_tbl_real = EorU_tbl_real[which(EorU_tbl_real$num_A < EorU_tbl_real$num_E),]

#write.table(Singlet_EorU_tbl_real, file = paste0("Table_", prefix, sample_name, "_Singlet_EorU_tbl_real_", run_number, sep = ""), quote = FALSE, sep = "\t")




###this is nothing
#Singlet_EorU_tbl_real[which(Singlet_EorU_tbl_real$num_edited > 10),]

###ok, so from here we'd like to:
#concatenate different tables
#rearrange columns in different tables so they match one another.
################
plate_id = "p3"

plateX_singlet_tbl = as.data.frame(read.table("/Users/anna/Documents/UW_GS/Shendure_Lab/Projects/Many_targets_ideas/CROPt/190722_Tree29_other7plates/TARG_plates2thru8/analysis/target_plate3/rerun_after_concatenation_fix/Table_190624_Tree29_plate3_Singlet_EorU_tbl_real_1"))
plateX_cigar_tbl = as.data.frame(read.table("/Users/anna/Documents/UW_GS/Shendure_Lab/Projects/Many_targets_ideas/CROPt/190722_Tree29_other7plates/TARG_plates2thru8/analysis/target_plate3/rerun_after_concatenation_fix/Table_190624_Tree29_plate3_Cigar_tbl_real_1"))
rownames(plateX_singlet_tbl) <- paste0(plate_id, rownames(plateX_singlet_tbl), sep = "")
rownames(plateX_cigar_tbl) <- paste0(plate_id, rownames(plateX_cigar_tbl), sep = "")

#for first plate
singlet_df = plateX_singlet_tbl
cigar_df = plateX_cigar_tbl

#for subsequent plates
singlet_df = rbind(singlet_df, plateX_singlet_tbl)
cigar_df = rbind(cigar_df, plateX_cigar_tbl)

#test_tbl_1 = cigar_df
#test_tbl_3 = plateX_cigar_tbl

#test_combined = rbind(test_tbl_1, test_tbl_3)
#test_combined[which(rownames(test_combined) == "p3A10_AACCTCAAGA"),]
#test_tbl_3[which(rownames(test_tbl_3) == "p3A10_AACCTCAAGA"),]

Singlet_EorU_tbl_real = singlet_df
Cigar_tbl_real = cigar_df

#setwd("/Users/anna/Documents/UW_GS/Shendure_Lab/Projects/Many_targets_ideas/CROPt/190722_Tree29_other7plates/TARG_plates2thru8/analysis/combined_p1_p3_take2")
#prefix = 
#write.table(Singlet_EorU_tbl_real, file = paste0("Table_", prefix, sample_name, "_Singlet_EorU_tbl_real_", run_number, sep = ""), quote = FALSE, sep = "\t")


#write.table(Singlet_EorU_tbl_real, file = "")

################

#only include cells w/ at least X edits for clustering
min_number_of_edits = 1
greater_than_edits = min_number_of_edits-1

Edited_singlet = Singlet_EorU_tbl_real[which(Singlet_EorU_tbl_real$num_edited > greater_than_edits),]
Edited_Cigar_tbl_real = Cigar_tbl_real[which(rownames(Cigar_tbl_real) %in% rownames(Edited_singlet)),]

library(plyr)
col_sums = as.data.frame(ldply(Edited_Cigar_tbl_real, function(c) sum(c=="X")))
col_sums$fraction_X = col_sums$V1/nrow(Edited_Cigar_tbl_real)
hist(col_sums$fraction_X, breaks = 50)
barplot(col_sums$fraction_X)
##check to make sure non between .5 & 1 -- if some are, change the .5 to something higher)
col_sums$fraction_X

only_less_than_50 = col_sums[which(col_sums$fraction_X < .8),]

only_colnames = colnames(Edited_Cigar_tbl_real)[colnames(Edited_Cigar_tbl_real) %in% only_less_than_50$.id]

Exclude_high_X_tbl = Edited_Cigar_tbl_real[,which(colnames(Edited_Cigar_tbl_real) %in% only_colnames)]

Edited_Cigar_tbl_real = Exclude_high_X_tbl

#Edited_EorU_tbl_real = EorU_tbl_real[which(EorU_tbl_real$num_edited > 0),]
#Edited_Cigar_tbl_real = Cigar_tbl_real[which(rownames(Cigar_tbl_real) %in% rownames(Edited_EorU_tbl_real)),]

#full_cigar_tbl = Cigar_tbl_real

temp_cigar_tbl_readl = Cigar_tbl_real
Cigar_tbl_real =temp_cigar_tbl_readl

Cigar_tbl_real = Edited_Cigar_tbl_real

##REAL TEMP SECTION
#temp_tbl = read.table("/Users/anna/Documents/UW_GS/Shendure_Lab/Projects/Many_targets_ideas/CROPt/190722_Tree29_other7plates/TARG_plates2thru8/analysis/all_plates/test_folder_for_tree_calculations/Plates_p1_p3_cigars_real_cells", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
#temp_blah = Cigar_tbl_real[,match(colnames(temp_tbl),colnames(Cigar_tbl_real))]
#Cigar_tbl_real = temp_blah
#write.table(Cigar_tbl_real, file = paste0("/Users/anna/Documents/UW_GS/Shendure_Lab/Projects/Many_targets_ideas/CROPt/190722_Tree29_other7plates/TARG_plates2thru8/analysis/combined_p1_p3_analysis/Cigar_tbl_for_corr_input_greater_than_", greater_than_edits,"_edited"), col.names = FALSE, quote = FALSE, sep = "\t")

#write.table(cbind(colnames(Cigar_tbl_real)), file = paste0("/Users/anna/Documents/UW_GS/Shendure_Lab/Projects/Many_targets_ideas/CROPt/190722_Tree29_other7plates/TARG_plates2thru8/analysis/combined_p1_p3_analysis/Cigar_tbl_for_corr_input_greater_than_", greater_than_edits,"_edited_COLNAMES"), col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
#write.table(rownames(Cigar_tbl_real), file = paste0("/Users/anna/Documents/UW_GS/Shendure_Lab/Projects/Many_targets_ideas/CROPt/190722_Tree29_other7plates/TARG_plates2thru8/analysis/combined_p1_p3_analysis/Cigar_tbl_for_corr_input_greater_than_", greater_than_edits,"_edited_ROWNAMES"), col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

#setwd("/Users/anna/Documents/UW_GS/Shendure_Lab/Projects/Many_targets_ideas/CROPt/190722_Tree29_other7plates/TARG_plates2thru8/analysis/combined_p1_p3_take2")
write.table(Cigar_tbl_real, file = paste0("Cigar_tbl_for_corr_input_greater_than_", greater_than_edits,"_edited"), col.names = FALSE, quote = FALSE, sep = "\t")
write.table(Cigar_tbl_real, file = paste0("Cigar_tbl_for_corr_input_greater_than_", greater_than_edits,"_edited_wCOLNAMES"), col.names = TRUE, quote = FALSE, sep = "\t")
write.table(cbind(colnames(Cigar_tbl_real)), file = paste0("Cigar_tbl_for_corr_input_greater_than_", greater_than_edits,"_edited_COLNAMES"), col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
write.table(rownames(Cigar_tbl_real), file = paste0("Cigar_tbl_for_corr_input_greater_than_", greater_than_edits,"_edited_ROWNAMES"), col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

##reorder columns based on past stuff;
tbl_w_cols_in_order = read.table("/Users/anna/Documents/UW_GS/Shendure_Lab/Projects/Many_targets_ideas/CROPt/191120_first_targ_atac/analysis/CIGAR_ordered_columns_w_names", header = TRUE, sep = "\t")
Cigar_tbl_real_ordered_cols = Cigar_tbl_real[,match(colnames(tbl_w_cols_in_order), colnames(Cigar_tbl_real))]
write.table(Cigar_tbl_real_ordered_cols, file = paste0("Cigar_tbl_for_corr_input_greater_than_", greater_than_edits,"_edited_ORDEREDwCOLNAMES"), col.names = TRUE, quote = FALSE, sep = "\t")

#### HERE RUN CPP SCRIPT 
##Maybe it's actually this one?? /Users/anna/Documents/UW_GS/Shendure_Lab/Projects/Many_targets_ideas/CROPt/CROPt_pipeline_files/Scripts/190702_distance_matrix_bash_wrapper_script.txt)
# ** cd into output folder!
# check that parameters are what you want them to be!
# /Users/anna/Documents/UW_GS/Shendure_Lab/Projects/Many_targets_ideas/CROPt/CROPt_pipeline_files/Scripts/1900820_distance_matrix_bash_wrapper_script_2.txt Cigar_tbl_for_corr_input_greater_than_0_edited

Run_num = 1 ##Have this match cpp script
#Output_dir = "/Users/anna/Documents/UW_GS/Shendure_Lab/Projects/Many_targets_ideas/CROPt/191126_targATAC_full_1/analysis/plate1/"
#setwd("/Users/anna/Documents/UW_GS/Shendure_Lab/Projects/Many_targets_ideas/CROPt/190722_Tree29_other7plates/TARG_plates2thru8/analysisNewPipeline/plate1")

normalized_dist_mat = read.table("Non_Normalized_Distance_Mat_1", header = FALSE, sep = "\t")
rownames(normalized_dist_mat) = normalized_dist_mat$V1
normalized_dist_mat = normalized_dist_mat[,c(2:ncol(normalized_dist_mat))]
colnames(normalized_dist_mat) = rownames(normalized_dist_mat)

#plot = pheatmap(normalized_dist_mat, fontsize = 1, silent = TRUE)
plot = pheatmap(normalized_dist_mat, clustering_method = "ward.D", fontsize = 1, silent = TRUE)
row_names_in_order = as.data.frame(rownames(normalized_dist_mat[plot$tree_row[["order"]],]))
colnames(row_names_in_order) = c("cells_in_order")
rownames(row_names_in_order) = row_names_in_order$cells_in_order

#pdf(paste0(Output_dir,"Heatmap_NON_normalized_Run_num", Run_num), width = 10, height = 10)
pdf(paste0("Heatmap_NON_normalized_Run_num", Run_num, "_wardD"), width = 10, height = 10)
plot
dev.off()

#plot2 = pheatmap(normalized_dist_mat, clustering_method = "ward.D", fontsize = 1)
#pdf(paste0(Output_dir,"Heatmap_NON_normalized_Run_num", Run_num, "_wardD"), width = 10, height = 10)
#plot2
#dev.off()

##### All this may be the subpar version of the stuff below
#EU_to_reorder = Del_or_Ins_EorU_tbl_real_ordered

#EU_reordered = EU_to_reorder[match(rownames(row_names_in_order), rownames(EU_to_reorder)),]
#EU_reordered_mat = as.matrix(EU_reordered)
#EU_reordered_mat_melt = melt(EU_reordered_mat)
#EU_reordered_mat_melt$value = as.factor(EU_reordered_mat_melt$value)
#EU_reordered_mat_melt$width = 1
#EU_reordered_mat_melt$height = .8

#plot_EU_reordered = ggplot(EU_reordered_mat_melt, aes(x = Var2, y = Var1)) +
#  geom_tile(aes(fill=value, width = EU_reordered_mat_melt$width, height = EU_reordered_mat_melt$height))
#cols = c("X" = "black", "D" = "RED", "U" = "yellow", "A" = "blue", "I" = "violet", "N" = "orange")

#plot_EU_reordered = plot_EU_reordered + scale_fill_manual(values = cols) + theme_minimal() + theme(legend.position="right") + theme(axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x=element_text(angle=90,hjust=1), axis.text.y=element_text(size = 6))

#pdf(paste0(Output_dir,"EU_reordered_", Run_num, ".pdf", sep = ""), width = 10, height = nrow(EU_reordered)/70*5)
#plot_EU_reordered
#dev.off()


####

cut_tree_1 = cutree(plot$tree_row, k = 1)
cut_tree_df = as.data.frame(cut_tree_1)

max_num_groups = nrow(normalized_dist_mat)

for(i in 2:max_num_groups)
{
  #col_name = paste0("cut_tree_", as.character(i))
  temp =  cutree(plot$tree_row, k = i)
  cut_tree_df = cbind(cut_tree_df, temp)
}
colnames(cut_tree_df) = c(1:max_num_groups)

write.table(cut_tree_df, file = paste0("cells_split_into_groups_wardD_Run_Num_", Run_num), quote = FALSE, sep = "\t")

#Cigar_tbl_real = as.data.frame(read.table("/Users/anna/Documents/UW_GS/Shendure_Lab/Projects/Many_targets_ideas/CROPt/190722_Tree29_other7plates/TARG_plates2thru8/analysis/combined_p1_p3_analysis/Cigar_tbl_for_corr_input_greater_than_0_edited", row.names = 1, header = FALSE))
#Cigar_tbl_header = read.table("/Users/anna/Documents/UW_GS/Shendure_Lab/Projects/Many_targets_ideas/CROPt/190722_Tree29_other7plates/TARG_plates2thru8/analysis/combined_p1_p3_analysis/Cigar_tbl_for_corr_input_greater_than_0_edited_COLNAMES", header = FALSE)
#colnames(Cigar_tbl_real) = Cigar_tbl_header$V1



cut_tree = cutree(plot$tree_row, k = 15)
#cut_tree_5 = cutree(plot$tree_row, k = 5)
#cut_tree_10 = cutree(plot$tree_row, k = 10)

#cut_tree_df = as.data.frame(cut_tree)
#cut_tree_df$cut_tree_5 = cut_tree_5
#colnames(cut_tree_df) = c("cut_tree_15", "cut_tree_5", "cut_tree_10")

plot(plot$tree_row)

plot$tree_row$merge

###for next section to work, rewrite cut_tree_df to just include one value


###TEST -- this is great!! makes a plot of edits indifferent colors
Cigar_tbl_real = as.data.frame(read.table("/Users/anna/Documents/UW_GS/Shendure_Lab/Projects/Many_targets_ideas/CROPt/190722_Tree29_other7plates/TARG_plates2thru8/analysis/combined_p1_p3_analysis/Cigar_tbl_for_corr_input_greater_than_0_edited", row.names = 1, header = FALSE))
Cigar_tbl_header = read.table("/Users/anna/Documents/UW_GS/Shendure_Lab/Projects/Many_targets_ideas/CROPt/190722_Tree29_other7plates/TARG_plates2thru8/analysis/combined_p1_p3_analysis/Cigar_tbl_for_corr_input_greater_than_0_edited_COLNAMES", header = FALSE)
colnames(Cigar_tbl_real) = Cigar_tbl_header$V1

EU_to_reorder = Cigar_tbl_real
EU_to_reorder$grey_col = "Blank"
#this part adds a last column -- you don't necessarily need to do this if you don't need it (skip to 601)
K = 20
cut_tree= cutree(plot$tree_row, k = K)
cut_tree_df = as.data.frame(cut_tree)
test = merge(EU_to_reorder, cut_tree_df["cut_tree"], by="row.names",all.x=TRUE)
rownames(test) = test$Row.names
test = test[,2:ncol(test)]
EU_to_reorder = test

EU_reordered = EU_to_reorder[match(rownames(row_names_in_order), rownames(EU_to_reorder)),]
EU_reordered_mat = as.matrix(EU_reordered)
EU_reordered_mat_melt = melt(EU_reordered_mat)
EU_reordered_mat_melt$value = as.factor(EU_reordered_mat_melt$value)
EU_reordered_mat_melt$width = 1
EU_reordered_mat_melt$height = .8

plot_EU_reordered = ggplot(EU_reordered_mat_melt, aes(x = Var2, y = Var1)) +
  geom_tile(aes(fill=value, width = EU_reordered_mat_melt$width, height = EU_reordered_mat_melt$height))
#cols = c("X" = "black", "D" = "RED", "U" = "yellow", "A" = "blue", "I" = "violet", "N" = "orange")

unique_cols = length(unique(EU_reordered_mat_melt$value))

set.seed(1)
cols = rainbow(unique_cols, s=.6, v=.9)[sample(1:unique_cols,unique_cols)]
new_tbl = as.data.frame(cbind(unique(as.character(EU_reordered_mat_melt$value))))
new_tbl$colors = cols
new_tbl[which(new_tbl$V1 == "70M:"),]$colors = "#FFFFFF"
new_tbl[which(new_tbl$V1 == "X"),]$colors = "black"
new_tbl[which(new_tbl$V1 == "AMB"),]$colors = "red"
new_tbl[which(new_tbl$V1 == "Blank"),]$colors = "gray"
#new_tbl$edit_color = paste("\"", as.character(new_tbl$V1),"\"=\"", as.character(new_tbl$colors), "\"", sep = "")
cols = as.character(new_tbl$colors)
names(cols) = as.character(new_tbl$V1)
#cols = c(as.character(as.factor(new_tbl$edit_color)))
#cols = c(test_list)

#test_list = paste0(new_tbl$edit_color, collapse=", ")

#names(new_tbl$V1) = levels(new_tbl$colors)

#plot_EU_reordered = plot_EU_reordered + scale_fill_manual(values = cols) + theme_minimal() + theme(legend.position="right") + theme(axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x=element_text(angle=90,hjust=1), axis.text.y=element_text(size = 6))
#plot_EU_reordered = plot_EU_reordered + scale_color_manual(name = new_tbl$V1, values = new_tbl$colors) + theme_minimal() + theme(legend.position="none") + theme(axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x=element_text(angle=90,hjust=1), axis.text.y=element_text(size = 6))
plot_EU_reordered = plot_EU_reordered + scale_fill_manual(values = cols) + theme_minimal() + theme(legend.position="none") + theme(axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x=element_text(angle=90,hjust=1), axis.text.y=element_text(size = 6))


#pdf(paste0(Output_dir,"Cigar_TEST_reordered_", Run_num, ".pdf", sep = ""), width = 10, height = nrow(EU_reordered)/70*5)
pdf(paste0(Output_dir,"Cigar_TEST_reordered_", Run_num, "WardD_K,", K, ".pdf", sep = ""), width = 10, height = nrow(EU_reordered)/70*5)
plot_EU_reordered
dev.off()

####Just one value using unique_cols from above

cut_tree= cutree(plot$tree_row, k = K)
cut_tree_df = as.data.frame(cut_tree)
cut_tree_df$cut_tree = as.numeric(as.character(cut_tree_df$cut_tree))

Special_output_dir = "/Users/anna/Documents/UW_GS/Shendure_Lab/Projects/Many_targets_ideas/CROPt/190722_Tree29_other7plates/TARG_plates2thru8/analysis/combined_p1_p3_take2/K20/"
for (i in 1:max(cut_tree_df$cut_tree))
{
  #i = 3
  temp_cut_tree_df = cut_tree_df
  temp_cut_tree_df$cell_names = rownames(temp_cut_tree_df)
  temp_cut_tree_df = temp_cut_tree_df[which(cut_tree_df$cut_tree == i),]
  temp_row_names_in_order = row_names_in_order[which(row_names_in_order$cells_in_order %in% rownames(temp_cut_tree_df)),]
  EU_to_reorder = Cigar_tbl_real[which(rownames(Cigar_tbl_real) %in% rownames(temp_cut_tree_df)),]
  EU_reordered = EU_to_reorder[match(temp_row_names_in_order, rownames(EU_to_reorder)),]
  print(head(EU_reordered))
  EU_reordered_mat = as.matrix(EU_reordered)
  EU_reordered_mat_melt = melt(EU_reordered_mat)
  EU_reordered_mat_melt$value = as.factor(EU_reordered_mat_melt$value)
  #print(nrow(EU_reordered_mat_melt))
  EU_reordered_mat_melt$width = 1
  EU_reordered_mat_melt$height = .8
  #unique_cols = length(unique(EU_reordered_mat_melt$value))
  plot_EU_reordered = ggplot(EU_reordered_mat_melt, aes(x = Var2, y = Var1)) +
    geom_tile(aes(fill=value, width = EU_reordered_mat_melt$width, height = EU_reordered_mat_melt$height))
  plot_EU_reordered = plot_EU_reordered + scale_fill_manual(values = cols) + theme_minimal() + theme(legend.position="none") + theme(axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x=element_text(angle=90,hjust=1), axis.text.y=element_text(size = 6))
  #plot_EU_reordered
  #pdf(paste0(Special_output_dir,"Cigar_TEST_reordered_", Run_num, "_K", max(cut_tree_df$cut_tree), "_cluster", i, ".pdf", sep = ""), width = 10, height = nrow(EU_reordered)/70*5)
  pdf(paste0(Special_output_dir,"Cigar_TEST_reordered_", Run_num, "_K", max(cut_tree_df$cut_tree), "_cluster", i, ".pdf", sep = ""), width = 10, height = nrow(EU_reordered)/70*5)
  print(plot_EU_reordered)
  dev.off()
}


###Make a new, smalled dataset from EU_reordered
subset_of_cells = c("p1B9_GATGCTACGA", "p1F10_ATGAACGCGC", "p3E12_GTCGACGGAA", "p1B10_TTCAAGAATC", "p3F7_TCCAGCAAT", "p3A2_TGCGCCTGGT", )
subset_EU_reordered = EU_reordered[which(),]





##### Getting clusters out https://stackoverflow.com/questions/6518133/clustering-list-for-hclust-function
#cut_tree = cutree(plot$tree_row, k = 25)
#cut_tree_df = as.data.frame(cut_tree)


#####





####This is all not necessary anymore


library(dplyr)
Cigar_tbl_real = sample_n(Cigar_tbl_real, 300)

#generate new correlation matrix
New_corr_matrix = as.data.frame(matrix(numeric(0),nrow(Cigar_tbl_real),nrow(Cigar_tbl_real)))
rownames(New_corr_matrix) = rownames(Cigar_tbl_real)
colnames(New_corr_matrix) = rownames(Cigar_tbl_real)

weight_same_edit = 5
weight_both_70M = 1
weight_one_X = 0
weight_at_least_one_AMB = 0
weight_diff_edits = -5
weight_one_70_one_edited = -2



for (i in 1:nrow(Cigar_tbl_real))
{
  print (i)
  for (j in 1:nrow(Cigar_tbl_real))
  {
    if (j>=i)
    {
      corr_weight = 0
      for (k in 1:ncol(Cigar_tbl_real))
      {
        if (Cigar_tbl_real[i,k] == "70M:" && Cigar_tbl_real[j,k] == "70M:")
        {
          corr_weight = corr_weight + weight_both_70M
        }
        else if (Cigar_tbl_real[i,k] == "X" || Cigar_tbl_real[j,k] == "X")
        {
          corr_weight = corr_weight + weight_one_X
        }
        else if (Cigar_tbl_real[i,k] == "AMB" || Cigar_tbl_real[j,k] == "AMB")
        {
          corr_weight = corr_weight + weight_at_least_one_AMB
        }
        else if (Cigar_tbl_real[i,k] == Cigar_tbl_real[j,k])
        {
          corr_weight = corr_weight + weight_same_edit
        }
        else if (Cigar_tbl_real[i,k] == "70M:" || Cigar_tbl_real[j,k] == "70M:")
        {
          corr_weight = corr_weight + weight_one_70_one_edited
        }
        else if (Cigar_tbl_real[i,k] != Cigar_tbl_real[j,k])
        {
          #print (paste0(Cigar_tbl_real[i,k],"---",Cigar_tbl_real[j,k]))
          corr_weight = corr_weight + weight_diff_edits
        }
      }
    }
    else
    {
      corr_weight = New_corr_matrix[j,i]
    }
    if(i == j)
    {
      #New_corr_matrix[i,j] = 0
      New_corr_matrix[i,j] = corr_weight
    }
    else
    {
      New_corr_matrix[i,j] = corr_weight
    }
  }
}

Non_neg_New_corr_matrix = New_corr_matrix + 77
Log_corr_mat = log(Non_neg_New_corr_matrix, base = 2)
log_plot = pheatmap(Log_corr_mat)

pdf(paste0(prefix, sample_name, "PHEATMAP_TEST_only_Edited_2_DIAG_ZEROS_600_cells_LOG_1.pdf", sep = ""), width = 100, height = 100)
log_plot
dev.off()

plot = pheatmap(New_corr_matrix)
row_names_in_order = as.data.frame(rownames(New_corr_matrix[plot$tree_row[["order"]],]))
write.table(row_names_in_order, file = "/Users/anna/Documents/UW_GS/Shendure_Lab/Projects/Many_targets_ideas/CROPt/190624_Tree29/analysis/300_DIAGNONZERO_1_tbl", col.names = FALSE, quote = FALSE)


pdf(paste0(prefix, sample_name, "PHEATMAP_TEST_only_Edited_2_DIAG_NONZEROS_300_cells_1.pdf", sep = ""), width = 50, height = 50)
plot
dev.off()


test_EU_tbl = EorU_tbl_real_ordered[which(rownames(EorU_tbl_real_ordered) %in% as.character(row_names_in_order$`rownames(New_corr_matrix[plot$tree_row[["order"]], ])`)),]
test_EU_tbl_reord = test_EU_tbl[match(rownames(test_EU_tbl),as.character(row_names_in_order$`rownames(New_corr_matrix[plot$tree_row[["order"]], ])`)),]



## THIS IS COPIES AND PASTED FROM ABOVE
EorU_tbl_real_ordered = as.data.frame(EorU_tbl_real_ordered_by_edits[,c(1:(ncol(EorU_tbl_real_ordered_by_edits) - 5))])
EorU_tbl_real_ordered_mat = as.matrix(EorU_tbl_real_ordered)
EorU_tbl_real_ordered_mat_melt = melt(EorU_tbl_real_ordered_mat)
EorU_tbl_real_ordered_mat_melt$value = as.factor(EorU_tbl_real_ordered_mat_melt$value)
EorU_tbl_real_ordered_mat_melt$width = 1
EorU_tbl_real_ordered_mat_melt$height = .8

plot_EorU_tbl_ordered = ggplot(EorU_tbl_real_ordered_mat_melt, aes(x = Var2, y = Var1)) +
  geom_tile(aes(fill=value, width = EorU_tbl_real_ordered_mat_melt$width, height = EorU_tbl_real_ordered_mat_melt$height))
cols = c("X" = "black", "E" = "RED", "U" = "yellow", "A" = "blue")

plot_EorU_tbl_ordered = plot_EorU_tbl_ordered + scale_fill_manual(values = cols) + theme_minimal() + theme(legend.position="right") + theme(axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x=element_text(angle=90,hjust=1), axis.text.y=element_text(size = 6))

pdf(paste0(prefix, sample_name, "_EorU_matrix_OREDERD_by_EA_plot_", run_number, ".pdf", sep = ""), width = 10, height = nrow(EorU_tbl_real)/70*5)
plot_EorU_tbl_ordered
dev.off()
##






#### ALL THIS IS JUST TESTING STUFF



Corr_tbl_real = Corr_tbl[which(rownames(Corr_tbl) %in% RNA_tbl_real$V1),which(colnames(Corr_tbl) %in% RNA_tbl_real$V1)]

m_test_melt = melt(as.matrix(Corr_tbl_real))
m_test_tbl_melt = m_test_melt
#m_test_tbl_melt$value = as.factor(m_test_tbl_melt$value)
m_test_tbl_melt$width = 1
m_test_tbl_melt$height = 1

ggplot(Corr_tbl) + geom_tile

plot_tbl = ggplot(m_test_melt, aes(x = Var2, y = Var1)) +
  geom_tile(aes(fill=value, width = m_test_melt$width, height = m_test_melt$height))

plot_tbl = ggplot(m_test_tbl_melt, aes(x = Var2, y = Var1)) +
  geom_tile(aes(fill=value, width = m_test_tbl_melt$width, height = m_test_tbl_melt$height), colour = "white") + scale_fill_gradient(low = "white", high = "steelblue")

pdf(paste0(prefix, sample_name, "_corr_plot_TEST", sep = ""), width = 20, height = 20)
plot_tbl
#bc_abun_plot = grid.arrange(plot_tbl, abunplot, widths = c(3,1), ncol=2)
dev.off()

most_corr = m_test_tbl_melt[order(-m_test_tbl_melt$value),]
most_corr_not_diag = most_corr[which(most_corr$Var1 != most_corr$Var2),]
most_corr_not_diag_ordered = most_corr_not_diag[order(-most_corr_not_diag$value ,most_corr_not_diag$Var1),]

Cigar_tbl[which(rownames(Cigar_tbl) %in% c("F2_CTCTTCAAGC")),]

Cigar_tbl[which(rownames(Cigar_tbl) %in% c("B5_GGCTTGCCAA","A4_CCGGCGGCGA", "F11_CATGAGAACT", "E6_CCGGCGGCGA", "D1_AAGTAATATT")),]
Cigar_tbl[which(rownames(Cigar_tbl) %in% c("F8_CATGAGAACT","H12_TCCAGCAATA", "C4_GGCTTGCCAA", "E8_TCCAGCAATA")),]




#heirarchical clustering test: THIS DOESN'T WORK WELL....
Cigar_tbl_real = Cigar_tbl[which(rownames(Cigar_tbl) %in% RNA_tbl_real$V1),]
install.packages("stats")

for_dist_mat = Corr_tbl_real

Normalized_corr_tbl = Corr_tbl_real
for(i in 1:nrow(Normalized_corr_tbl))
{
  Normalized_corr_tbl[i,] = Normalized_corr_tbl[i,]-Normalized_corr_tbl[i,i]
}

Normalized_corr_tbl = -1*Normalized_corr_tbl
Normalized_corr_tbl_NA = Normalized_corr_tbl


for (i in 1:nrow(Normalized_corr_tbl_NA))
{
  for(j in 1:nrow(Normalized_corr_tbl_NA))
  {
    if (j>i)
    {
      Normalized_corr_tbl_NA[i,j] = "NA"
    }
  }
}

Normalized_corr_tbl_NA = as.dist(as.matrix(Normalized_corr_tbl_NA, diag = TRUE))


plot(hclust(Normalized_corr_tbl_NA,  method = "complete", members = NULL))
plot(hclust(Normalized_corr_tbl_NA))
pdf(paste0(prefix, sample_name, "HEIRARCHICAL_TEST", sep = ""), width = 100, height = 20)
plot(hclust(Normalized_corr_tbl_NA))
dev.off()

Cigar_tbl_real[which(rownames(Cigar_tbl_real) %in% c("A6_AAGTTACCTA","F12_TGCAGCCTAC","A2_GGAGTAAGCC")),]

install.packages("pheatmap")
library(pheatmap)

plot = pheatmap(Corr_tbl_real)

pdf(paste0(prefix, sample_name, "PHEATMAP_TEST.pdf", sep = ""), width = 50, height = 50)
plot
dev.off()

Normalized_corr_tbl = Corr_tbl_real
for(i in 1:nrow(Normalized_corr_tbl))
{
  Normalized_corr_tbl[i,] = Normalized_corr_tbl[i,]-Normalized_corr_tbl[i,i]
}

plot_norm = pheatmap(Normalized_corr_tbl)
pdf(paste0(prefix, sample_name, "PHEATMAP_TEST_NORMALIZED.pdf", sep = ""), width = 50, height = 50)
plot_norm
dev.off()

