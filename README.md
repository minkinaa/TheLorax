# TheLorax

Data & Code Availabiliy

## Raw & Processed Data:

RNA, on GEO:
- Raw expression data (fastq), separated by plates. Each plate contains 96 Read1-Read2 paired files, each associated with a single well (i.e. index #2 during combinatorial indexing). Read2 contains a UMI & index #1. An index #1 list is available on github.
- Processed cell by gene matrix, including raw counts

ATAC, on GEO:
- Raw atac data (fastq), with sequences renamed by single cell identifiers (e.g. p2G1_TTAACGCCGT-R1.3 is associated with cell G1_TTAACGCCGT from plate2; these identified match those found in processed lineage data).
- Processed cell by gene matrix (includes 5kb upstream of promoter), including raw counts per interval

Lineage_wRNA, on GEO:
Fastqs associated with lineage targets captured alongside RNA, separated by plate & well/index #2

Lineage_wATAC, on GEO:
Fastqs associated with lineage targets captured alongside ATAC, separated by plate & well/index #2 (plate 8 does not exist)

Processed Lineage Data (github):
lineage_profiles_wATAC.txt
lineage_profiles_wRNA.txt
tree_file_LinATAC
tree_file_LinRNA


## Code:

RNA processing to generate cell by gene files implemented in wrapper script 190807_sciRNA_wrapper_ALL.txt. Internal code & input files have been uploaded to the folder "RNA_Processing" but may contain hard-coded file paths.
