#' Kaplan-Meier plots
#' 
#' This Rscript demonstrates joining BigQuery tables to compare survival rates conditional on somatic mutations. 
#' 
#required packages
require(bigrquery) || install.packages('bigrquery')
library(survival)
source('ggsurv.R')

# specify cloud project name
cloudProject = "isb-cgc"

# BigQuery dataset you want to work with 
dataset = "tcga_201510_alpha"

#list tables in the dataset
bigrquery::list_tables(cloudProject,dataset)

#we will use Somatic_Mutation_calls and Clinical_data for this example
somatic_mutation_table = paste("[",cloudProject,":",dataset,".Somatic_Mutation_calls",']',sep="")
clinical_table = paste("[",cloudProject,":",dataset,".Clinical_data",']',sep="")

#1. First, let's see which protein change is most common across all tumors. Let's limit the result to 20 most frequent protein changes. 
# (GROUP BY ParticipantBarcode and Protein_Change to eliminate duplicate entries)
sqlQuery = paste("SELECT Protein_Change,count(*) as N",
                 "FROM (SELECT ParticipantBarcode,Protein_Change",
                 "FROM", somatic_mutation_table,
                 "GROUP BY ParticipantBarcode,Protein_Change)",
                 "WHERE Protein_Change <> 'null'",
                 "GROUP BY Protein_Change",
                 "ORDER BY N DESC",
                 "LIMIT 20")
result = query_exec(sqlQuery,project = cloudProject)
#save result as a data frame
resultDF = data.frame(result)
#view results
resultDF


#2. Now, for the most frequent protein change, let's count the number of samples per tumor type with this protein change
#Note: The example protein changes used here are unique to their genes, but generally a protein change can be associated with
#more than one genes

top_protein_change = resultDF[1,1]

sqlQuery = paste("SELECT Study, count(*) as n ",
                 "FROM (SELECT ParticipantBarcode, Study ", 
                  "FROM ",somatic_mutation_table,
                  " WHERE Protein_Change='",top_protein_change,"' ",   
                  "GROUP BY ParticipantBarcode,Protein_Change, Study) ",
                  "GROUP BY Study ORDER BY n DESC",sep="")
result = query_exec(sqlQuery,project = cloudProject)
resultDF = data.frame(result)
resultDF

#let's look at only top 10 studies for survival comparison.
top_ten_studies = resultDF$Study[1:10]


#3. let's see how survival trends look across different tumor types with the most common protein change
sqlQuery = paste("SELECT clin.ParticipantBarcode,vital_status, clin.Study, days_to_last_known_alive ",
                 "FROM ",clinical_table," as clin JOIN (SELECT ParticipantBarcode, Study ", 
                 "FROM ",somatic_mutation_table,
                 " WHERE Protein_Change='",top_protein_change,"' ",   
                 "GROUP BY ParticipantBarcode,Protein_Change, Study) as mut ", 
                 "ON clin.ParticipantBarcode=mut.ParticipantBarcode",sep="")
result = query_exec(sqlQuery,project = cloudProject)
resultDF = data.frame(result)
head(resultDF)

#convert vital status to numeric type
resultDF$vital_status[resultDF$vital_status=="Alive"] = 0
resultDF$vital_status[resultDF$vital_status=="Dead"] = 1
class(resultDF$vital_status) = "numeric"

#convert Study type to factor
resultDF$clin_Study = factor(resultDF$clin_Study)

#normalize days_to_last_known_alive within each tumor type to mitigate inherent survival differences across tumor types
for (tumor_type in levels(resultDF$clin_Study))
{
  survival_this_tumor = resultDF$days_to_last_known_alive[resultDF$clin_Study==tumor_type]
  
  min_last_known_alive = min(survival_this_tumor,na.rm = TRUE)
  max_last_known_alive = max(survival_this_tumor,na.rm = TRUE)
  
  resultDF$days_to_last_known_alive[resultDF$clin_Study==tumor_type] = (resultDF$days_to_last_known_alive[resultDF$clin_Study==tumor_type] - min_last_known_alive)/(max_last_known_alive - min_last_known_alive)
  
}

#fit survival curve to data
tcga.surv = survfit(Surv(resultDF$days_to_last_known_alive,resultDF$vital_status)~resultDF$clin_Study,data = resultDF)
#plot survival curve
ggsurv(tcga.surv,main = paste("Protein Change = ",top_protein_change))


#4. Let's only plot top_ten_studies (if there are more than 10 curves in the previous plot)
resultDF = resultDF[resultDF$clin_Study %in% top_ten_studies,]
#fit survival curve to data
tcga.surv = survfit(Surv(resultDF$days_to_last_known_alive,resultDF$vital_status)~resultDF$clin_Study,data = resultDF)
#plot survival curve
ggsurv(tcga.surv,main = paste("Protein Change = ",top_protein_change))


#5. let's add a baseline survival curve to this plot. For this, we'll pull all the samples from top_ten_studies that do not have 
#the top_protein_change
#Make the following simple changes to the query above
#* LEFT JOIN clinical_table to somatic_mutation_table, so that the query result contains all the samples from the clinical data
# regardles of whether they have a mutation or not.
#* Add Protein_Change to list of output columns. It's value will be null for non-matching samples.
#* Add a WHERE clause to get samples only from top_ten_studies

sqlQuery = paste("SELECT clin.ParticipantBarcode,vital_status, clin.Study, days_to_last_known_alive, Protein_Change ",
                 "FROM ",clinical_table," as clin LEFT JOIN (SELECT ParticipantBarcode, Protein_Change, Study ", 
                 "FROM ",somatic_mutation_table,
                 " WHERE Protein_Change='",top_protein_change,"' ",   
                 "GROUP BY ParticipantBarcode,Protein_Change, Study) as mut ", 
                 "ON clin.ParticipantBarcode=mut.ParticipantBarcode ",
                 "WHERE clin.Study in (",paste(shQuote(top_ten_studies),collapse = ","),")",sep="")
result = query_exec(sqlQuery,project = cloudProject)
resultDF = data.frame(result)
head(resultDF)

#convert vital status to numeric type
resultDF$vital_status[resultDF$vital_status=="Alive"] = 0
resultDF$vital_status[resultDF$vital_status=="Dead"] = 1
class(resultDF$vital_status) = "numeric"

#make all controls to be in 'All_Controls' study type
resultDF$clin_Study[is.na(resultDF$Protein_Change)] = "All_Controls"

#convert Study type to factor
resultDF$clin_Study = factor(resultDF$clin_Study)

#normalize days_to_last_known_alive within each tumor type to mitigate inherent survival differences across tumor types
for (tumor_type in unique(resultDF$clin_Study))
{
  survival_this_tumor = resultDF$days_to_last_known_alive[resultDF$clin_Study==tumor_type]
  
  min_last_known_alive = min(survival_this_tumor,na.rm = TRUE)
  max_last_known_alive = max(survival_this_tumor,na.rm = TRUE)
  
  resultDF$days_to_last_known_alive[resultDF$clin_Study==tumor_type] = (resultDF$days_to_last_known_alive[resultDF$clin_Study==tumor_type] - min_last_known_alive)/(max_last_known_alive - min_last_known_alive)
  
}

#fit survival curve to data
tcga.surv = survfit(Surv(resultDF$days_to_last_known_alive,resultDF$vital_status)~resultDF$clin_Study,data = resultDF)
#plot survival curve
ggsurv(tcga.surv,main = paste("Protein Change = ",top_protein_change))


#6. Let's look at the survival curve for the tumor type where the top_protein_change is most common.
#For this, we need one simple change to the previous sql query.
#* edit the WHERE claus to return only top_ten_studies[1] records
# Note: Also try top_protein_change='p.V600E' and top_study='THCA' 

sqlQuery = paste("SELECT clin.ParticipantBarcode,vital_status, clin.Study, days_to_last_known_alive, Protein_Change ",
                 "FROM ",clinical_table," as clin LEFT JOIN (SELECT ParticipantBarcode, Protein_Change, Study ", 
                 "FROM ",somatic_mutation_table,
                 " WHERE Protein_Change='",top_protein_change,"'",   
                 " GROUP BY ParticipantBarcode,Protein_Change, Study) as mut", 
                 " ON clin.ParticipantBarcode=mut.ParticipantBarcode",
                 " WHERE clin.Study='",top_ten_studies[1],"'",sep="")

result = query_exec(sqlQuery,project = cloudProject)
resultDF = data.frame(result)
head(resultDF)

#convert vital status to numeric type
resultDF$vital_status[resultDF$vital_status=="Alive"] = 0
resultDF$vital_status[resultDF$vital_status=="Dead"] = 1
class(resultDF$vital_status) = "numeric"

#Convert Protein_Change column to a TRUE/FALSE column so we can compare survival trends between cases and controls
resultDF$Protein_Change = as.logical(!is.na(resultDF$Protein_Change))

#fit survival curve to data
tcga.surv = survfit(Surv(resultDF$days_to_last_known_alive,resultDF$vital_status)~resultDF$Protein_Change,data = resultDF)
#plot survival curve
ggsurv(tcga.surv,main = paste("Protein Change = ",top_protein_change,", Study = ",top_ten_studies[1]))

