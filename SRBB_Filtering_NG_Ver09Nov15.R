SRBFF_Filtering <- function(input_count_table_file,input_gene_names_file,output_directory,min_total_each_treatment,min_samples_with_min_reads,min_reads_per_sample_group)
{
	count_table <- read.csv(input_count_table_file,header=T)
	print("Input matrix")
	print(head(count_table))
	
	print("**************************")
	print(c("Total genes of input count table: ",nrow(count_table)))
	print("**************************")
	
    gene_names <- read.csv(input_gene_names_file,header=T)
    
    row.names(count_table) <- gene_names[,1]
    
	flag_treat_reps<-"n"
	
	treatment_groups <- NA
	
	while(is.na(treatment_groups)==TRUE)
	{
		treatment_groups <- as.numeric(readline("Please enter the total treatment groups:"))  
	}
	
	treat_reps_matrix <- matrix(nrow=2,ncol=treatment_groups,data=0)
    row.names(treat_reps_matrix) <- c("Treatments_Bio_Reps","Total_Reads")
	
	while (flag_treat_reps=="n")
	{
				
		flag_equ_reps <- as.character(readline("Is the number of biological replicates equal for all treatment groupt? y=yes, n=no   "))
		
		if ((flag_equ_reps=="y") || (flag_equ_reps=="yes"))
		{
			reps_treatment_groups <- as.numeric(readline("Please enter the total biological replicates for each treatment groups: "))
			
			total_treat_reps <- treatment_groups * reps_treatment_groups
			
			print(c("The count table includes total of ",total_treat_reps,"samples (data columns),based on", treatment_groups, "treatment groups and ", reps_treatment_groups, "replicates per treatment"),quote=F)			

			treat_reps_matrix [1,] <- reps_treatment_groups
        }
		
		else if (flag_equ_reps=="n"|| (flag_equ_reps=="no"))
		
		{
			print("Please enter the number of biological replicates for each treatment group.")
			

			for(c_treat in 1:treatment_groups)
			{
				print(c("Treatment group number: ", c_treat),quote=F)
				treat_reps_matrix[1,c_treat]<- as.numeric(readline("Please enter the total biological replicates for treatment group: "))
			}			
		}
		
		
		if(anyNA(treat_reps_matrix))
		{
			print("The biological replicate numbers for all treatment groups should be specified! Let's try again!")
		}
		else
		{
			flag_treat_reps<-"y"
		}
	}
			
			

	
    Output_matrix <- matrix(data = 1:ncol(count_table),nrow=1,ncol=ncol(count_table))
	colnames(Output_matrix) <- colnames(count_table)
	
	Tossed_matrix <- matrix(data = 1:ncol(count_table),nrow=1,ncol=ncol(count_table))
	colnames(Tossed_matrix) <- colnames(count_table)
	
	flag_stageI <- 0
	flag_stageII <- 0
	
	
    for(c_rows in 1:nrow(count_table))
	{
		#Please use the following commented line(s) to track the progress
		#print(c("Working on gene number: ",c_rows))
		
		begin_col <- 1
		for(c_bio_reps in 1:ncol(treat_reps_matrix))
		{
			for(c_cols in begin_col:(begin_col+treat_reps_matrix[1,c_bio_reps]-1))
				{
					treat_reps_matrix[2,c_bio_reps] <- treat_reps_matrix[2,c_bio_reps] + count_table[c_rows,(c_cols)]
				}
			begin_col <- begin_col + treat_reps_matrix[1,c_bio_reps]
			
		}
	
		#Please use the following commented line(s) to track the progress
		#print("Total reads for each treatment group (for all its biological replicates):")
		#print(treat_reps_matrix)

				
        if (any((as.numeric(treat_reps_matrix[2,]) >= as.numeric(min_total_each_treatment))))
		{
			#Please use the following commented line(s) to track the progress
			#print("Stage I: Passed!")	
			
			flag_stageI <- flag_stageI + 1
			Output_matrix <- rbind(Output_matrix,as.matrix(count_table[c_rows,]))
			

		}
		else 
		{
			flag_min_reads_per_sample_group <- 0
			for(c_sums in 1:(ncol(treat_reps_matrix)))
			{
				if(as.numeric(treat_reps_matrix[2,c_sums]) >= as.numeric(min_reads_per_sample_group))
				{
                    flag_min_reads_per_sample_group <- flag_min_reads_per_sample_group + 1
				}
			}
            
	
			if(as.numeric(flag_min_reads_per_sample_group) >= as.numeric(min_samples_with_min_reads))
			{
				#Please use the following commented line(s) to track the progress
				#print("Stage II: Passed!")
				
				flag_stageII <- flag_stageII + 1		
				Output_matrix <- rbind(Output_matrix,count_table[c_rows,])
			}
			else
			{
				#Please use the following commented line(s) to track the progress
				#print("Stage II: Did NOT Pass!")
				
                Tossed_matrix <- rbind(Tossed_matrix,count_table[c_rows,])
				
			}
            
            if(anyNA(Tossed_matrix))
            {
                print("Tossed_matrix error")
                exit()
            }
            if(anyNA(Output_matrix))
            {
                print("Output_matrix error")
                exit()
            }
 
		}
	
        treat_reps_matrix [2,] <- 0
        flag_min_reads_per_sample_group <- 0
        
		#Please use the following commented line(s) to track the progress
		#print("*****************************************************************************************")
	}
	
	print("**************************")
    print(c("Total processed genes: ",nrow(count_table)))
	print(c("Total genes of kept genes: ",nrow(Output_matrix) - 1))
	print(c("Total genes of kept genes, in Stage I of filtering: ",flag_stageI))
	print(c("Total genes of kept genes, in Stage II of filtering: ",flag_stageII))
	print("---------")
    print(c("Total genes of filtered out genes: ",nrow(Tossed_matrix) - 1))
	print("**************************")

	write.csv(Output_matrix[2:nrow(Output_matrix),],paste(output_directory,paste("Kept_Genes_",min_total_each_treatment,"_",min_samples_with_min_reads,"_",min_reads_per_sample_group,"_",Sys.time(),".csv",sep=""),sep=""))

	write.csv(Tossed_matrix[2:nrow(Tossed_matrix),],paste(output_directory,paste("Tossed_Genes_",min_total_each_treatment,"_",min_samples_with_min_reads,"_",min_reads_per_sample_group,"_",Sys.time(),".csv",sep=""),sep=""))
	
}