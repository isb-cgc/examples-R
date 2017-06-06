# Working with barcode lists

As part of working with the ISB-CGC web app, you will have created cohorts, represented by lists
of barcodes. This short tutorial shows how to retrieve your cohorts, query for sample details,
and compose a BigQuery, all from within the R environment.

The isb-cgc project has a collection of web services called 'endpoints', which accept
and return information. The endpoints allow the user to interact with the isb-cgc system programmatically,
or for a client application to do so on behalf of the user.

The ISBCGCExamples package contains a number of "wrapper" functions, that make calling
the endpoints possible from the R environment. These are just examples, and much
more is possible. See the documentation to get more ideas:

http://isb-cancer-genomics-cloud.readthedocs.io/en/latest/sections/progapi/Programmatic-API.html#isb-cgc-api

To get started, we load up R, and import the ISBCGCExamples library.


```r
library(ISBCGCExamples)
library(bigrquery)
```

### Creating a token

The first step is creating a token. This token contains your authentication
status, and lets the service know about what information is available to you.


```r
my_token <- isb_init()
```

Calling the isb_init function is going to open a browser window that lets you
authenticate with Google.

### Listing cohorts

To get a listing of the previously created cohorts, we can use the list_cohorts
function that takes a token, and returns a list with items including 'count',
'items', 'kind', and 'etag'. The count shows the number of saved cohorts and
the items contains information about the cohorts.


```r
my_cohorts <- list_cohorts(my_token)
names(my_cohorts)
```

```
## [1] "count" "items" "kind"  "etag"
```

If there is a count of four, then the my_cohorts$items list will have length four.
Each item has a name. We can use the lapply function (list-apply) to get the
names out of each list element.


```r
lapply(my_cohorts$items, function(x) x$name)
```

```
##[[1]]
##[1] "All TCGA Data"

##[[2]]
##[1] "HPV_cohort"

##[[[3]]
##[[1] "GBM_Adult_TP"

##[[[4]]
##[[1] "BRCA_Adult_TP"

##[[[5]]
##[[1] "HNSC_Adult_TP_and_NT"
```

Also, importantly, each element in 'items' has an 'id', which is used when we
do further queries on the cohort.


```r
lapply(my_cohorts$items, function(x) x$id)
```

```
##[[[1]]
##[[1] "1"

##[[[2]]
##[[1] "403"

##[[[3]]
##[[1] "847"

##[[[4]]
##[[1] "848"

##[[[5]]
##[[1] "849"
```

### Getting barcode lists from a cohorts

Now that we have the cohort IDs, we can collect the various barcodes contained
in the cohort. These include patient barcodes, sample barcodes, and platform
specific aliquot barcodes. To do this, we can use the barcodes_from_cohort function.


```r
my_cohort_id <- lapply(my_cohorts$items, function(x) x$id)[[4]]
my_barcodes <- cohort_barcodes(my_cohort_id, my_token)
names(my_barcodes)
```

```
## [1] "name"            "sample_count"    "permission"      "source_type"
## [5] "samples"         "comments"        "id"              "parent_id"
## [9] "patients"        "source_notes"    "filters"         "last_date_saved"
## [13] "email"           "patient_count"   "kind"            "etag"
```

The object returned from barcodes_from_cohort is again a list, this time with
elements 'cohort_id', 'sample_count', 'patient_count', 'patients', and 'samples'.
The patients and samples elements are also lists, but lists of patients or
sample barcodes.


```r
my_barcodes$cases[1:5]
```

```
##[[[1]]
##[[1] "TCGA-E2-A10C"

##[[[2]]
##[[1] "TCGA-B6-A0I8"

[[3]]
[1] "TCGA-C8-A278"

[[4]]
[1] "TCGA-C8-A1HI"

[[5]]
[1] "TCGA-AO-A0J8"
```

```r
my_barcodes$samples[1:5]
```

```
## [[1]]
## [1] "TCGA-2W-A8YY-01A"
##
## [[2]]
## [1] "TCGA-4J-AA1J-01A"
##
## [[3]]
## [1] "TCGA-4P-AA8J-01A"
##
## [[4]]
## [1] "TCGA-BA-4074-01A"
##
## [[5]]
## [1] "TCGA-BA-4075-01A"
```

### Getting details on a sample

Suppose you have a particular sample of interest. We can use the endpoints to
get details about the sample.


```r
my_sample_barcode <- my_barcodes$samples[[1]]
my_sample_details  <- sample_get(my_sample_barcode)

names(my_sample_details)
```

```
## [1] "data_details"       "patient"            "data_details_count"
## [4] "aliquots"           "biospecimen_data"   "kind"
## [7] "etag"
```

```r
my_sample_details$data_details_count
```

```
## [1] 9
```

```r
length(my_sample_details$data_details)
```

```
## [1] 9
```

Here the data_details_count list element tells us that there are 18 "data details"
which contain information like data platforms, file names, paths to cloud storage,
the TCGA data level.

For example, let's see what platforms are represented.


```r
lapply(my_sample_details$data_details, function(x) paste(x$platform, ", ", x$experimental_strategy, sep=""))
```

```
##[[1]]
##[[1] "Clinical,"

##[[2]]
##[[1] "Clinical,"

##[[3]]
##[[1] "Illumina GA,miRNA-Seq"

##[[4]]
##[[1] "Illumina GA,WXS"

##[[5]]
##[[1] "Illumina HiSeq,RNA-Seq"

##[[6]]
##[[1] "Illumina GA,miRNA-Seq"

##[[7]]
##[[1] "Illumina,miRNA-Seq"

##[[8]]
##[[1] "Illumina,WXS"

##[[9]]
##[[1] "Illumina,RNA-Seq"
```


### Moving from the endpoints API to BigQuery

If we have a list of barcodes, then we can programmatically construct a BigQuery.

First let's flatten the list of sample barcodes to a character vector.


```r
# get the character vector of samples
samples <- unlist(my_barcodes$samples)

# SQL strings need to be surrounded by single quotes
samples_with_quotes <- sapply(samples[1:10], function(x) paste("'",x,"'", sep=""))

# then the samples need to be surrounded by parenthesis.
query_samples <- paste("( ", paste(samples_with_quotes, collapse=","), " )")

bq <- paste("
SELECT
   project_short_name,
   AVG(LOG2(normalized_count+1)) as mean_log2_expression
FROM
   [isb-cgc:TCGA_hg19_data_v0.RNAseq_Gene_Expression_UNC_RSEM]
WHERE
   HGNC_gene_symbol = 'EGFR' AND
   sample_barcode IN ", query_samples,"
GROUP BY
   project_short_name"
)

results <- query_exec(bq, project)
results
```

```
##project_short_name mean_log2_expression
##1          TCGA-BRCA             6.668811
```

At this point, if this is your first query, a browser window will pop-up, and
you will need to authenticate for BigQuery.

Since we're working in R, we can take advantage of the great visualization libraries.

### Wrap up

To reiterate, using the endpoints API, we got the list of previously created cohorts, retrieved the
list of sample barcodes, examined data associated with the samples, and constructed
a BigQuery.

More information on the endpoints API and R package example functions can be found at:

http://isb-cancer-genomics-cloud.readthedocs.io/en/latest/sections/progapi/Programmatic-API.html

https://github.com/isb-cgc/examples-R/blob/master/inst/doc/Working_With_Barcode_Lists.md




```r
sessionInfo()
```

```
##R version 3.3.2 (2016-10-31)
##Platform: x86_64-apple-darwin13.4.0 (64-bit)
##Running under: OS X El Capitan 10.11.6
##
##locale:
##[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
##
##attached base packages:
##[1] stats     graphics  grDevices utils     datasets  methods   base
##
##other attached packages:
##[1] httr_1.2.1           bigrquery_0.3.0.9000 ISBCGCExamples_0.1.3
##
##loaded via a namespace (and not attached):
## [1] magrittr_1.5   R6_2.2.0       assertthat_0.1 tools_3.3.2    DBI_0.6
## [6] dplyr_0.5.0    curl_2.4       tibble_1.2     Rcpp_0.12.11   jsonlite_1.3
##[11] httpuv_1.3.3   openssl_0.9.6
```
