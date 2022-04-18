# TheLorax

Data & Code Availabiliy

# GEO:

RNA:
- Raw expression data (fastq), separated by plates. Each plate contains 96 Read1-Read2 paired files, each associated with a single well (i.e. index #2 during combinatorial indexing). Read2 contains a UMI & index #1. An index #1 list is available on github.
- Processed cell by gene matrix, including raw counts

ATAC:
- Raw atac data (fastq), with sequences renamed by single cell identifiers (e.g. p2G1_TTAACGCCGT-R1.3 is associated with cell G1_TTAACGCCGT from plate2; these identified match those found in processed lineage data).
- Processed cell by gene matrix (includes 5kb upstream of promoter), including raw counts per interval

Lineage_wRNA:
Fastqs associated with lineage targets captured alongside RNA, separated by plate & well/index #2

Lineage_wATAC:
Fastqs associated with lineage targets captured alongside ATAC, separated by plate & well/index #2 (plate 8 does not exist)


## Code:
