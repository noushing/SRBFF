This document describes the parameters needed for running general-purpose filtering method Sample-Replicate based Feature Filtering (SRBFF). The code is developed by Noushin Ghaffari at Texas A&M University. The SRBFF is provided as a R function and here is an example of command to call the SRBFF function:

SRBFF_Filtering(Count_Table_Reads.csv, Count_Table_Gene_Names.csv, Output_Directory_Path, min_read_for_treatments_in_StageI, treatments_with_min_reads_StageII, min_read_per_treatments_in_StageII)

The input parameters of the SRBBFF function are:

1) Count_Table_Reads.csv is a comma-delimited text file (CSV), where each column has the treatments and their biological replicates, and each row represents the features (e.g. genes or exons), and each cell holds the reads mapped to each sample for each feature. First line is the header line that has the name of treatments. Here is an example:

Treatment1_Replicate1, Treatment1_Replicate2, Treatment2_Replicate1, Treatment2_Replicate2
10,15,1,1
3,0,4,9
...

The count table should have at least two treatments with at least 1 biological replicates.

2) Count_Table_Gene_Names.csv is a comma-delimited text file that has a header row, and the rest of the rows hold the name of the features. The number of rows in the Count_Table_Reads.csv file and this file should match exactly. If the name of the features is not known, one can simply supply a text file that has dummy values for gene names. Here is an example of Count_Table_Gene_Names.csv file:

Gene_names
G1
G2
G3
...

3) Output_Directory_Path is the absolute path to the directory to save the output files, which includes a CSV file for the kept genes (Filtered_Genes_*.csv) and genes that have been discarded (Tossed_Genes_*.csv).

4)  The value "min_read_for_treatments_in_StageI" defines the total number of reads that is required in stage I of SRBFF (please refer to the SRBFF description for more details). This value is the summation of the reads mapped to all biological replicates of at least one treatment to keep the feature in stage I.

5  and 6) The "treatments_with_min_reads_StageII" is the minimum conditions (the treatment with all its replicates) that are required to have the "min_read_per_treatments_in_StageII" to pass stage II of SRBFF.

As an example, let's assume "min_read_for_treatments_in_StageI, treatments_with_min_reads_StageII, min_read_per_treatments_in_StageII" are 20, 2, 15 respectively. At stage I all the genes that have at least 20 reads for one treatment (all biological reps included) are kept. Genes that didn't pass stage I, considered in stage II. If at least 2 treatments have total of 15 mapped reads (all biological reps included), then, the gene passes stage II and is kept. All other genes are discarded.



 
