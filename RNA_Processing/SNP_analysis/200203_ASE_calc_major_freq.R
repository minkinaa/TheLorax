#This is the first script for determining major & minor allele from data

library(ggplot2)

return_letter = function(x){
  if (x[1] == x[5]){
    return ("A")
  }
  else if(x[2] == x[5]){
    return ("C")
  }
  else if(x[3] == x[5]){
    return ("G")
  }
  else if(x[4] == x[5]){
    return ("T")
  }
}

return_major_count = function(x){
  if (x[1] == x[5]){
    return (x[1])
  }
  else if(x[2] == x[5]){
    return (x[2])
  }
  else if(x[3] == x[5]){
    return (x[3])
  }
  else if(x[4] == x[5]){
    return (x[4])
  }
}

second_highest = function(x){
  temp_vals = vector()
  for(i in 1:4){
    if (x[i] != x[5])
    {
      temp_vals = c(temp_vals,x[i])
    }
  }
  return(max(temp_vals))
}
  
return_count_from_letter = function(x)
{
  if(x[5] == "A")
  {
    return(x[1])
  }
  else if(x[5] == "C"){
    return(x[2])
  }
  else if(x[5] == "G"){
    return(x[3])
  }
  else if(x[5] == "T"){
    return(x[4])
  }
}

output_dir = "/Users/anna/Documents/UW_GS/Shendure_Lab/Projects/Many_targets_ideas/CROPt/190722_Tree29_other7plates/TARG_plates2thru8/200623_NewPipelineAnalysis/210326_ModTree42_RNA_Analysis/211109_RNA_chr_copy_number_plots/"

name_seq = c(as.character(1:22),"X")

#for(i in 1:length(name_seq)){
  #chr = paste0("chr", name_seq[i])
chr = "chr14"

chrom_ref_file = read.table("/Users/anna/Documents/UW_GS/Shendure_Lab/Projects/Many_targets_ideas/CROPt/CROPt_pipeline_files/Commonly_sourced_files/chrom_length_midpoint_etc", sep = "\t", stringsAsFactors = FALSE, header = T)
chrom_ref_file$chr_num = as.character(chrom_ref_file$chr_num)
chrom_ref_file[which(chrom_ref_file$chr_num == "23"),]$chr_num = "X"

all_cells_ASE_tbl = read.table(paste0("/Users/anna/Documents/UW_GS/Shendure_Lab/Projects/Many_targets_ideas/CROPt/190722_Tree29_other7plates/TARG_plates2thru8/200623_NewPipelineAnalysis/210326_ModTree42_RNA_Analysis/RNA_ASE_Analysis/AllPlates_",chr,"_grp1_42_base_counts_per_position_wSNPinfo"), header = TRUE, stringsAsFactors = FALSE)

all_cells_ASE_tbl$sum = all_cells_ASE_tbl$A + all_cells_ASE_tbl$C + all_cells_ASE_tbl$T + all_cells_ASE_tbl$G
all_cells_ASE_tbl$max = apply(all_cells_ASE_tbl[,c(3:6)], 1, max)
all_cells_ASE_tbl$frac_max = all_cells_ASE_tbl$max/all_cells_ASE_tbl$sum
all_cells_ASE_tbl_HET = all_cells_ASE_tbl[which(all_cells_ASE_tbl$frac_max < .85 & all_cells_ASE_tbl$SNP != "none"),]

hist(all_cells_ASE_tbl[which(all_cells_ASE_tbl$SNP != "none"),]$frac_max, breaks = 50)
hist(all_cells_ASE_tbl_HET$frac_max, breaks = 30)

all_cells_ASE_tbl_HET$exp_major = apply(all_cells_ASE_tbl_HET[,c(3:6,11)], 1, return_letter)
all_cells_ASE_tbl_HET$exp_major_count = apply(all_cells_ASE_tbl_HET[,c(3:6,11)], 1, return_major_count)

all_cells_ASE_tbl_HET$exp_minor_count = apply(all_cells_ASE_tbl_HET[,c(3:6,11)], 1, second_highest)
all_cells_ASE_tbl_HET$exp_minor = apply(all_cells_ASE_tbl_HET[,c(3:6,15)], 1, return_letter)

plot(all_cells_ASE_tbl_HET$pos, all_cells_ASE_tbl_HET$frac_max, pch = 20)

chr_length = as.numeric(chrom_ref_file[which(paste0("chr",chrom_ref_file$chr_num) == chr),]$chr_length)
centomere_pos = (as.numeric(chrom_ref_file[which(paste0("chr",chrom_ref_file$chr_num) == chr),]$cent_end) - as.numeric(chrom_ref_file[which(paste0("chr",chrom_ref_file$chr_num) == chr),]$cent_start))/2 + as.numeric(chrom_ref_file[which(paste0("chr",chrom_ref_file$chr_num) == chr),]$cent_start)
breakpoint = 0

hist(all_cells_ASE_tbl_HET$frac_max, breaks = 30, main = chr)

all_cells_ASE_tbl_HET_left = all_cells_ASE_tbl_HET[which(all_cells_ASE_tbl_HET$pos < breakpoint),]
all_cells_ASE_tbl_HET_right = all_cells_ASE_tbl_HET[which(all_cells_ASE_tbl_HET$pos > breakpoint),]

sum(all_cells_ASE_tbl_HET_left$exp_minor_count)/sum(all_cells_ASE_tbl_HET_left$exp_major_count)
sum(all_cells_ASE_tbl_HET_right$exp_minor_count)/sum(all_cells_ASE_tbl_HET_right$exp_major_count)

bin_size = 5000000
bin_num = round(chr_length/bin_size)

minor_ov_major = vector()
midpoint_vect = vector()
major_over_total_vect = vector()
total_counts_major = vector()
total_counts_all = vector()
total_num_values_per_bin = vector()
for(i in 1:(bin_num)){
  temp_start = i*bin_size - bin_size
  temp_end = temp_start+bin_size
  temp_midpoint = temp_start + (bin_size/2)
  temp_tbl = all_cells_ASE_tbl_HET[which(all_cells_ASE_tbl_HET$pos > temp_start & all_cells_ASE_tbl_HET$pos < temp_end),]
  print(nrow(temp_tbl))
  temp_ratio = sum(temp_tbl$exp_minor_count)/sum(temp_tbl$exp_major_count)
  temp_major_over_total = sum(temp_tbl$exp_major_count)/(sum(temp_tbl$exp_minor_count) + sum(temp_tbl$exp_major_count))
  minor_ov_major = c(minor_ov_major, temp_ratio)
  midpoint_vect = c(midpoint_vect, temp_midpoint)
  major_over_total_vect = c(major_over_total_vect, temp_major_over_total)
  total_counts_major = c(total_counts_major, sum(temp_tbl$exp_major_count))
  total_counts_all = c(total_counts_all, sum(temp_tbl$exp_major_count) + sum(temp_tbl$exp_minor_count))
  total_num_values_per_bin = c(total_num_values_per_bin, nrow(temp_tbl))
}

bin_count_tbl = cbind(midpoint_vect,minor_ov_major)
plot(bin_count_tbl)

total_counts_all_LARGE = as.data.frame(cbind(c(1:length(total_counts_all)), total_counts_all))

bin_count_tbl_major_ov_total = cbind(midpoint_vect, major_over_total_vect, total_num_values_per_bin)
plot(bin_count_tbl_major_ov_total, main = paste0(chr, ", bin_size=", bin_size), ylim = c(.5,1),  xlab = paste0(chr, " pos, bin_size=", bin_size), ylab = "Major allele frequency", pch = 20)
abline(v = breakpoint)
abline(h = .5)
abline(h = .666)
abline(h = .75)
abline(h = 3/5)
abline(h = .55)


bin_count_tbl_major_ov_total = as.data.frame(bin_count_tbl_major_ov_total)

bin_count_tbl_major_ov_total$log2_size = log(bin_count_tbl_major_ov_total$total_num_values_per_bin, base = 2)

#### add ggplot here
scatt_plot_w_scale = ggplot(as.data.frame(bin_count_tbl_major_ov_total)) +
  theme_bw() + ylim(.5,1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_hline(yintercept = 1/2, linetype="dashed", color = "gray") +
  geom_hline(yintercept = 2/3, linetype="dashed", color = "gray") +
  geom_hline(yintercept = 3/4, linetype="dashed", color = "gray") +
  geom_hline(yintercept = 1, linetype="dashed", color = "gray") +
  geom_vline(xintercept = centomere_pos, linetype="solid", color = "yellow") +
  geom_point(aes(x = midpoint_vect, y = major_over_total_vect, col = log2_size), size = 2) +
  geom_point(aes(x = midpoint_vect, y = major_over_total_vect), shape = 1, size = 2, colour = "black") +
  theme(axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank()) +
  #scale_color_gradient2(low="white", high="black", limits = c(0, 10), oob = scales::squish, name = "")
  scale_colour_gradient2(low="darkblue", mid = "white", high="darkred", limits = c(0, 10), midpoint = 5, oob = scales::squish, name = "")
 
scatt_plot_NO_scale = ggplot(as.data.frame(bin_count_tbl_major_ov_total)) +
  theme_bw() + ylim(.5,1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_hline(yintercept = 1/2, linetype="dashed", color = "gray") +
  geom_hline(yintercept = 2/3, linetype="dashed", color = "gray") +
  geom_hline(yintercept = 3/4, linetype="dashed", color = "gray") +
  geom_hline(yintercept = 1, linetype="dashed", color = "gray") +
  geom_vline(xintercept = centomere_pos, linetype="solid", color = "yellow") +
  geom_point(aes(x = midpoint_vect, y = major_over_total_vect, col = log2_size), size = 2) +
  geom_point(aes(x = midpoint_vect, y = major_over_total_vect), shape = 1, size = 2, colour = "black") +
  theme(axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank()) +
  #scale_color_gradient2(low="white", high="black", limits = c(0, 10), oob = scales::squish, name = "")
  scale_colour_gradient2(low="darkblue", mid = "white", high="darkred", limits = c(0, 10), midpoint = 5, oob = scales::squish, name = "") +
  theme(legend.position = "none")

scatt_plot_NO_scale_wCHR_name = ggplot(as.data.frame(bin_count_tbl_major_ov_total)) +
  theme_bw() + ylim(.5,1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_hline(yintercept = 1/2, linetype="dashed", color = "gray") +
  geom_hline(yintercept = 2/3, linetype="dashed", color = "gray") +
  geom_hline(yintercept = 3/4, linetype="dashed", color = "gray") +
  geom_hline(yintercept = 1, linetype="dashed", color = "gray") +
  geom_vline(xintercept = centomere_pos, linetype="solid", color = "yellow") +
  geom_point(aes(x = midpoint_vect, y = major_over_total_vect, col = log2_size), size = 2) +
  geom_point(aes(x = midpoint_vect, y = major_over_total_vect), shape = 1, size = 2, colour = "black") +
  theme(axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank()) +
  #scale_color_gradient2(low="white", high="black", limits = c(0, 10), oob = scales::squish, name = "")
  scale_colour_gradient2(low="darkblue", mid = "white", high="darkred", limits = c(0, 10), midpoint = 5, oob = scales::squish, name = "") +
  theme(legend.position = "none") + ggtitle(chr)


#pdf(paste0(output_dir,paste0("wName_DotPlotMajorAlleleFreq_",chr,"binsize",format(bin_size, scientific=F),".pdf"), sep = ""), width = 3.5, height = 2)
#print(scatt_plot_NO_scale_wCHR_name)
#dev.off()

#}
  
#pdf(paste0(output_dir,paste0("DotPlotMajorAlleleFreq_",chr,"binsize",format(bin_size, scientific=F),"_wScale.pdf"), sep = ""), width = 3.5, height = 2)
#scatt_plot_w_scale
#dev.off()

pdf(paste0(output_dir,paste0("DotPlotMajorAlleleFreq_",chr,"binsize",format(bin_size, scientific=F),".pdf"), sep = ""), width = 3.5, height = 2)
scatt_plot_NO_scale
dev.off()

## A test here:
temp_tbl$total_count = temp_tbl$A + temp_tbl$C + temp_tbl$G + temp_tbl$T
het_temp_tbl = temp_tbl[which(temp_tbl$major_count != temp_tbl$total_count),]
nrow(het_temp_tbl[which(het_temp_tbl$minor_count > 0),])
nrow(het_temp_tbl[which(het_temp_tbl$minor_count == 0),])



###Let's make some nice code for this part.

##Step 1: make a new global table from all cells:
chr = "chr3"
group = "grp10:19"

all_cells_ASE_tbl = read.table(paste0("/Users/anna/Documents/UW_GS/Shendure_Lab/Projects/Many_targets_ideas/CROPt/190722_Tree29_other7plates/TARG_plates2thru8/200623_NewPipelineAnalysis/210326_ModTree42_RNA_Analysis/RNA_ASE_Analysis/AllPlates_",chr,"_grp1_42_base_counts_per_position_wSNPinfo"), header = TRUE, stringsAsFactors = FALSE)
grp1_ASE_tbl_name = "/Users/anna/Documents/UW_GS/Shendure_Lab/Projects/Many_targets_ideas/CROPt/190722_Tree29_other7plates/TARG_plates2thru8/200623_NewPipelineAnalysis/210326_ModTree42_RNA_Analysis/RNA_ASE_Analysis/AllPlates_chr3_grp10_19_base_counts_per_position"
grp1_ASE_tbl = read.table(grp1_ASE_tbl_name, header = TRUE, stringsAsFactors = FALSE)

chrom_ref_file = read.table("/Users/anna/Documents/UW_GS/Shendure_Lab/Projects/Many_targets_ideas/CROPt/CROPt_pipeline_files/Commonly_sourced_files/chrom_length_midpoint_etc", sep = "\t", stringsAsFactors = FALSE, header = T)
chrom_ref_file$chr_num = as.character(chrom_ref_file$chr_num)
chrom_ref_file[which(chrom_ref_file$chr_num == "23"),]$chr_num = "X"

chr_length = as.numeric(chrom_ref_file[which(paste0("chr",chrom_ref_file$chr_num) == chr),]$chr_length)
centomere_pos = (as.numeric(chrom_ref_file[which(paste0("chr",chrom_ref_file$chr_num) == chr),]$cent_end) - as.numeric(chrom_ref_file[which(paste0("chr",chrom_ref_file$chr_num) == chr),]$cent_start))/2 + as.numeric(chrom_ref_file[which(paste0("chr",chrom_ref_file$chr_num) == chr),]$cent_start)

bin_size = 5000000
bin_num = round(chr_length/bin_size)


all_cells_ASE_tbl$sum = all_cells_ASE_tbl$A + all_cells_ASE_tbl$C + all_cells_ASE_tbl$T + all_cells_ASE_tbl$G
all_cells_ASE_tbl$max = apply(all_cells_ASE_tbl[,c(3:6)], 1, max)
all_cells_ASE_tbl$frac_max = all_cells_ASE_tbl$max/all_cells_ASE_tbl$sum
all_cells_ASE_tbl_HET = all_cells_ASE_tbl[which(all_cells_ASE_tbl$frac_max < .85 & all_cells_ASE_tbl$SNP != "none"),]

#hist(all_cells_ASE_tbl[which(all_cells_ASE_tbl$SNP != "none"),]$frac_max, breaks = 50)
#hist(all_cells_ASE_tbl_HET$frac_max, breaks = 30)

all_cells_ASE_tbl_HET$exp_major = apply(all_cells_ASE_tbl_HET[,c(3:6,11)], 1, return_letter)
all_cells_ASE_tbl_HET$exp_major_count = apply(all_cells_ASE_tbl_HET[,c(3:6,11)], 1, return_major_count)

all_cells_ASE_tbl_HET$exp_minor_count = apply(all_cells_ASE_tbl_HET[,c(3:6,11)], 1, second_highest)
all_cells_ASE_tbl_HET$exp_minor = apply(all_cells_ASE_tbl_HET[,c(3:6,15)], 1, return_letter)

## Now read the group file
#grp1_ASE_tbl = read.table(grp1_ASE_tbl_name, header = TRUE, stringsAsFactors = FALSE)
grp1_ASE_tbl = grp1_ASE_tbl[which(grp1_ASE_tbl$pos %in% all_cells_ASE_tbl_HET$pos),]
temp_large_table_subset = all_cells_ASE_tbl_HET[which(paste0(all_cells_ASE_tbl_HET$chr, "_", all_cells_ASE_tbl_HET$pos) %in% paste0(grp1_ASE_tbl$chr, "_", grp1_ASE_tbl$pos)),]

grp1_ASE_tbl$exp_major = temp_large_table_subset$exp_major
grp1_ASE_tbl$exp_minor = temp_large_table_subset$exp_minor
grp1_ASE_tbl$major_count = apply(grp1_ASE_tbl[,c(3:6,7)], 1, return_count_from_letter)
grp1_ASE_tbl$minor_count = apply(grp1_ASE_tbl[,c(3:6,8)], 1, return_count_from_letter)
grp1_ASE_tbl$major_count = as.numeric(as.character(grp1_ASE_tbl$major_count))
grp1_ASE_tbl$minor_count = as.numeric(as.character(grp1_ASE_tbl$minor_count))
grp1_ASE_tbl$two_counts = 0
grp1_ASE_tbl[which(grp1_ASE_tbl$major_count > 0 & grp1_ASE_tbl$minor_count > 0),]$two_counts = 1
#plot(grp1_ASE_tbl$pos, grp1_ASE_tbl$two_counts, pch = 20)

minor_ov_major = vector()
midpoint_vect = vector()
major_over_total_vect = vector()
total_counts_major = vector()
total_counts_all = vector()
total_num_values_per_bin = vector()
for(i in 1:(bin_num)){
#for(i in 1:8){
  temp_start = i*bin_size - bin_size
  temp_end = temp_start+bin_size
  temp_midpoint = temp_start + (bin_size/2)
  temp_tbl =grp1_ASE_tbl[which(grp1_ASE_tbl$pos > temp_start & grp1_ASE_tbl$pos < temp_end),]
  print(nrow(temp_tbl))
  temp_ratio = sum(temp_tbl$minor_count)/sum(temp_tbl$major_count)
  temp_major_over_total = sum(temp_tbl$major_count)/(sum(temp_tbl$minor_count) + sum(temp_tbl$major_count))
  minor_ov_major = c(minor_ov_major, temp_ratio)
  midpoint_vect = c(midpoint_vect, temp_midpoint)
  major_over_total_vect = c(major_over_total_vect, temp_major_over_total)
  total_counts_major = c(total_counts_major, sum(temp_tbl$major_count))
  total_counts_all = c(total_counts_all, sum(temp_tbl$major_count) + sum(temp_tbl$minor_count))
  total_num_values_per_bin = c(total_num_values_per_bin, nrow(temp_tbl))
}

bin_count_tbl = cbind(midpoint_vect,minor_ov_major)
plot(bin_count_tbl)

bin_count_tbl_major_ov_total = cbind(midpoint_vect, major_over_total_vect, total_num_values_per_bin)
#total_counts_all_LARGE$small_group = total_counts_all
#total_counts_all_LARGE$percent = total_counts_all_LARGE$small_group/total_counts_all_LARGE$total_counts_all
#plot(total_counts_all_LARGE$V1, total_counts_all_LARGE$percent/median(total_counts_all_LARGE[which(!is.nan(total_counts_all_LARGE$percent)),]$percent), main = group)

bin_count_tbl_major_ov_total = as.data.frame(bin_count_tbl_major_ov_total)
plot(bin_count_tbl_major_ov_total, main = paste0(chr, ", bin_size=", bin_size, ", ", group), ylim = c(.5,1), xlab = paste0(chr, " pos, bin_size=", bin_size), ylab = "Major allele frequency", pch = 20)
abline(v = breakpoint)
#abline(v = 21593751)
#abline(v = 100881329)
#abline(v = 58500000)
abline(h = .666)
abline(h = .75)
abline(h = .5)

bin_count_tbl_major_ov_total = as.data.frame(bin_count_tbl_major_ov_total)

bin_count_tbl_major_ov_total$log2_size = log(bin_count_tbl_major_ov_total$total_num_values_per_bin, base = 2)

scatt_plot_NO_scale_wCHR_name = ggplot(as.data.frame(bin_count_tbl_major_ov_total)) +
  theme_bw() + ylim(0,1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_hline(yintercept = 1/2, linetype="dashed", color = "gray") +
  geom_hline(yintercept = 2/3, linetype="dashed", color = "gray") +
  geom_hline(yintercept = 3/4, linetype="dashed", color = "gray") +
  geom_hline(yintercept = 1, linetype="dashed", color = "gray") +
  geom_vline(xintercept = centomere_pos, linetype="solid", color = "yellow") +
  geom_point(aes(x = midpoint_vect, y = major_over_total_vect, col = log2_size), size = 2) +
  geom_point(aes(x = midpoint_vect, y = major_over_total_vect), shape = 1, size = 2, colour = "black") +
  theme(axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank()) +
  #scale_color_gradient2(low="white", high="black", limits = c(0, 10), oob = scales::squish, name = "")
  scale_colour_gradient2(low="darkblue", mid = "white", high="darkred", limits = c(0, 10), midpoint = 5, oob = scales::squish, name = "") +
  theme(legend.position = "none") + ggtitle(paste0(chr, ", ", group))

pdf(paste0(output_dir,paste0("DotPlotMajorAlleleFreq_",chr,"binsize",format(bin_size, scientific=F),"_",group,".pdf"), sep = ""), width = 3.5, height = 2)
scatt_plot_NO_scale_wCHR_name
dev.off()

scatt_plot_NO_scale_wCHR_name


###let's explore what is happening at pos 8: when major allele not dominant, is minor allele present (or noise)?
##temp_tbl generated in loop that ends at bin 8 (randomly chosen)
##this is chr6, where one of 4 copies is lost (we are saying) and everything becomes homozygous; so all should be major allele
temp_tbl$total_count = temp_tbl$A + temp_tbl$C + temp_tbl$G + temp_tbl$T
het_temp_tbl = temp_tbl[which(temp_tbl$major_count != temp_tbl$total_count),]
nrow(het_temp_tbl[which(het_temp_tbl$minor_count > 0),])
nrow(het_temp_tbl[which(het_temp_tbl$minor_count == 0),])


###ok, let's calculate noise vs. not....

library(ggplot2)

return_whether_major_base_expected = function(x){
  expected_bases = as.character(x[c(8,9)])
  if(x[12] %in% expected_bases){
    return("yes")
  } else {
    return("no")
  }
}

return_whether_minor_base_expected = function(x){
  expected_bases = as.character(x[c(8,9)])
  if(x[15] %in% expected_bases){
    return("yes")
  } else {
    return("no")
  }
}

name_seq = c(as.character(1:22),"X")

percent_major_vect = vector()
percent_minor_vect = vector()

for(i in 1:length(name_seq)){
  chr = paste0("chr",name_seq[i])
  print(chr)

  all_cells_ASE_tbl = read.table(paste0("/Users/anna/Documents/UW_GS/Shendure_Lab/Projects/Many_targets_ideas/CROPt/190722_Tree29_other7plates/TARG_plates2thru8/200623_NewPipelineAnalysis/210326_ModTree42_RNA_Analysis/RNA_ASE_Analysis/AllPlates_",chr,"_grp1_42_base_counts_per_position_wSNPinfo"), header = TRUE, stringsAsFactors = FALSE)
  just_snps = all_cells_ASE_tbl[which(all_cells_ASE_tbl$SNP != "none"),]
  just_snps_simple = just_snps[which(just_snps$major %in% c("A", "T", "C", "G") & just_snps$minor %in% c("A", "T", "C", "G")),]
  
  just_snps_simple$sum = just_snps_simple$A + just_snps_simple$C + just_snps_simple$T + just_snps_simple$G
  just_snps_simple$max = apply(just_snps_simple[,c(3:6)], 1, max)
  
  just_snps_simple$exp_major = apply(just_snps_simple[,c(3:6,11)], 1, return_letter)
  just_snps_simple$exp_major_count = apply(just_snps_simple[,c(3:6,11)], 1, return_major_count)
  
  just_snps_simple$exp_minor_count = apply(just_snps_simple[,c(3:6,11)], 1, second_highest)
  just_snps_simple$exp_minor = apply(just_snps_simple[,c(3:6,14)], 1, return_letter)
  
  just_snps_simple$major_expected = apply(just_snps_simple, 1, return_whether_major_base_expected)
  just_snps_simple$minor_expected = apply(just_snps_simple, 1, return_whether_minor_base_expected)
  
  percent_major_expected = nrow(just_snps_simple[which(just_snps_simple$major_expected == "yes"),])/nrow(just_snps_simple)*100
  percent_minor_expected = nrow(just_snps_simple[which(just_snps_simple$minor_expected == "yes"),])/nrow(just_snps_simple)*100
  
  percent_major_vect = c(percent_major_vect,percent_major_expected)
  percent_minor_vect = c(percent_minor_vect, percent_minor_expected)

}

all_chr_tbl = as.data.frame(cbind(name_seq,percent_major_vect,percent_minor_vect))

hist(as.numeric(all_chr_tbl$percent_minor_vect))
barplot(as.numeric(all_chr_tbl$percent_minor_vect))

all_chr_tbl$name_seq <- factor(all_chr_tbl$name_seq, levels=name_seq)

minor_plot = ggplot(all_chr_tbl, aes(x = name_seq, y = as.numeric(percent_minor_vect))) + 
  geom_bar(stat="identity", width = .8, fill = "skyblue") +
  theme_classic() +
  theme(axis.title.y=element_blank(), axis.title.x=element_blank())

major_plot = ggplot(all_chr_tbl, aes(x = name_seq, y = as.numeric(percent_major_vect))) + 
  geom_bar(stat="identity", width = .8, fill = "salmon") +
  theme_classic() +
  theme(axis.title.y=element_blank(), axis.title.x=element_blank())


output_dir = "/Users/anna/Documents/UW_GS/Shendure_Lab/Projects/Many_targets_ideas/CROPt/190722_Tree29_other7plates/TARG_plates2thru8/200623_NewPipelineAnalysis/210326_ModTree42_RNA_Analysis/211109_RNA_chr_copy_number_plots/"

pdf(paste0(output_dir,paste0("Percent_minor_expected.pdf"), sep = ""), width = 4.5, height = 2)
minor_plot
dev.off()

pdf(paste0(output_dir,paste0("Percent_major_expected.pdf"), sep = ""), width = 4.5, height = 2)
major_plot
dev.off()








