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



