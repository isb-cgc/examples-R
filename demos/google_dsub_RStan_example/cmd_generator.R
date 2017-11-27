
# writes out a task matrix for dsub.

files <- list.files("data")

#--env SAMPLE_ID --input DATA_FILE    --output OUTPUT_TABLE   --output OUTPUT_PLOT
tasks <- data.frame(check.names = FALSE)

bucket <- "gs://gibbs_bucket_nov162016/"

for (i in 1:length(files)) {
    print(i)
    sampleid <- i
    input_file <- paste(bucket,"data/",files[i],sep="")
    output_plot <- paste(bucket,"stan_plot",i,".png",sep="")
    output_table <- paste(bucket,"stan_table",i,".txt",sep="")
    tasks <- rbind(tasks, data.frame('--env SAMPLE_ID'=sampleid, '--input DATA_FILE'=input_file, '--output OUTPUT_TABLE'=output_table, '--output OUTPUT_PLOT'=output_plot, check.names = FALSE))
}

write.table(tasks, file="task_matrix.txt", row.names=F, col.names=T, sep='\t', quote=F)
