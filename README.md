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

### Cell IDs found in the processed data file are related to raw fastq files as follows:

sciRNA_lib: Each fastq pair file contains information from 25 cells (or less). The cell IDs are a combination of the plate/batch, PCR well index (index #2) during combinatorial indexing (in fastq file name), and index #1, the last 10 bases of the the read1 sequence. E.g. the raw data for the cell called "p1A10_ACTAATTGAG", can be found in the raw data file "p1A10-RNA_CCAGGTCTAC_S10", with those sequences with "ACTAATTGAG" as the last 10 bases in the R1 file corresponding to this cell. The first 8 bases of R1 are the UMI.

sciATAC_lib: Cell IDs corresponding to those found in the processed file can be found in the sequence name in the fastq (e.g. "@p2H1_TGAGCTACTT-R1.1" --> cell ID = "p2H1_TGAGCTACTT")

Lineage_wRNA & Lineage_wATAC: Cell IDs are related to fastq file names as described for "sciRNA_lib" except index #1 is located in positions 9-18 in the R1 fastq file. The first 8 bases of R1 are the UMI. Read structure is further describe in manuscript Materials & Methods.


## Code:

*Some scripts contain hard-coded file paths

RNA processing to generate cell by gene files implemented in wrapper script 190807_sciRNA_wrapper_ALL.txt. Internal code & input files have been uploaded to the folder "RNA_Processing." Downstream processing and visualization code located in internal folders. See manuscript for further details. 

Initial lineage data processing to generate cell lineage profiles is outlined in Lineage_Generate_Profiles_Sample_Commands.txt. Internal code and input files in folder "Lineage_Processing." See manuscript for further details.

Code and associated files for downstream lineage processing (splitting duplicated targets, correcting missing data, tree building & visualization) can be found in folder "Lineage_Tree_Building_etc." See manuscript for further details.




