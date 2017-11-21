

#gcloud alpha genomics pipelines run \
#--pipeline-file standocker-pipeline.yaml \
#--inputs INPUT_FILE=gs://gibbs_bucket_nov162016/data/data_file_1.csv \
#--inputs INPUT_SCRIPT=gs://gibbs_bucket_nov162016/logistic_regression_ref_man.R \
#--outputs OUTPUT_PLOT=gs://gibbs_bucket_nov162016/stan_output/stan_test_plot1.png \
#--outputs OUTPUT_FILE=gs://gibbs_bucket_nov162016/stan_output/stan_test_table1.txt \
#--logging gs://gibbs_bucket_nov162016/logs/


library(glue)

files <- read.table("file_list.txt")

cmdlist <- data.frame()

bucket <- "gs://your-google-bucket-name/"

pipeline_file <- "standocker-pipeline.yaml"
input_script <- "logistic_regression_ref_man.R"

for (i in 1:nrow(files)) {
    print(i)
     input_file <- files[i,1]
     output_plot <- paste("stan_test_plot",i,".png",sep="")
     output_table <- paste("stan_test_table",i,".txt",sep="")
     cmd <- glue("gcloud alpha genomics pipelines run --pipeline-file {pipeline_file}  --inputs INPUT_FILE={bucket}data/{input_file} --inputs INPUT_SCRIPT={bucket}{input_script} --outputs OUTPUT_PLOT={bucket}stan_output/{output_plot} --outputs OUTPUT_FILE={bucket}stan_output/{output_table} --logging {bucket}logs/")
     cmdlist <- rbind(cmdlist, data.frame(CMD=cmd))
}

write.table(cmdlist, file="cmds.txt", row.names=F, col.names=F, quote=F)
