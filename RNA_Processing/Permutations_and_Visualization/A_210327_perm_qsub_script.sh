num_iterations=$1
run_num=$2

Rscript \
/net/shendure/vol1/home/minkinaa/scripts/191203_NEW_CROPt_pipeline_files/201216_CLUSTER_smarter_exp_group_permutations_R_faster.R \
/net/shendure/vol10/projects/AM_Heterogeneity/nobackup/190722_Tree29_all_plates/210326_Tree29_42Grp_RNA_NewAnalysis/210326_1001_per_10_times_faster_script/ \
/net/shendure/vol1/home/minkinaa/reference_files/unique_gencode.v29.gene_ids_gene_names_PLUS.bed \
/net/shendure/vol10/projects/AM_Heterogeneity/nobackup/190722_Tree29_all_plates/210326_Tree29_42Grp_RNA_NewAnalysis/All_Cell_GroupNum_and_ConsCell_Lookup_NewestGroupsAdded \
/net/shendure/vol10/projects/AM_Heterogeneity/nobackup/190722_Tree29_all_plates/201125_Tree29_RNA_NewAnalysis/All_NonLog_Normalized_Scaled_gene_counts \
${num_iterations} \
${run_num} \
/net/shendure/vol10/projects/AM_Heterogeneity/nobackup/190722_Tree29_all_plates/210326_Tree29_42Grp_RNA_NewAnalysis/ReorderedCutDF_aka_cutdf

