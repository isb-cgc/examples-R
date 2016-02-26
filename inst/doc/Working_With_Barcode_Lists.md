# Working with lists of barcodes

Supposing we have gcloud installed and the authorization process is completed
with an active isb-cgc project ID.

First from the command line we need to download the scripts needed to
authorize with ISB and access file lists.

Documentation is found here: http://isb-cancer-genomics-cloud.readthedocs.org/en/latest/sections/progapi/Programmatic-API.html

### Saving time
If you just want info on a couple barcodes, then doing them one at a time is fine.
But, once you have a large set of barcodes, it's much faster to create a cohort
with that list of barcodes, and query on the whole cohort together.
*Coming soon*!

```
#<< ON the command line >>#

# the isb authorization scrips need the oauth2client library
[sudo] pip install --upgrade oauth2client

# download and run the isb authorization script.
wget https://raw.githubusercontent.com/isb-cgc/ISB-CGC-Webapp/master/scripts/isb_auth.py
python isb_auth.py --noauth_local_webserver

# then download the isb script to access file lists.
wget https://raw.githubusercontent.com/isb-cgc/ISB-CGC-Webapp/master/scripts/isb_curl.py
```

Now start up an R session.

```
# handy string processing library to parse the file list returns.
library(stringr)


parseCurl <- function(x) {
  # function takes the return string from isb_curl (x)
  # and returns a list of file paths.
  y <- str_trim(x)
  idx <- str_detect(pattern="^\"gs://", string=y)  # lines containing files
  z <- str_split(y[idx], "\"")
  unlist(lapply(z, function(zi) zi[2]))
}

callCurl <- function(barcode) {
  # here we build the command line command to query the web service with isb_curl.py
  print(paste0("working on: ", barcode))
  cmd <- paste0("python isb_curl.py https://api-dot-isb-cgc.appspot.com/_ah/api/cohort_api/v1/datafilenamekey_list_from_sample?sample_barcode=",barcode)
  system(cmd, intern=T)
}


# read in a list of barcodes
barcodes <- read.table("sample_barcode_list.txt", stringsAsFactors=F)
barcodes <- unique(barcodes$V1)

# we can get sample details on a single sample barcode
cmd <- paste0("python isb_curl.py https://api-dot-isb-cgc.appspot.com/_ah/api/cohort_api/v1/sample_details?sample_barcode=TCGA-44-7670-10A")
system(cmd, intern=T)

{
 "data_details": [
  {
   "Project": "TCGA",
   "Pipeline": "broad.mit.edu__snp_cnv",
   "SampleBarcode": "TCGA-44-7670-10A",
   "Repository": "DCC",
   "DatafileUploaded": "false",
   "Datatype": "Copy Number Results-SNP",
```
   etc ....

We can also get patient details using the patient barcode ...

```
cmd <- paste0("python isb_curl.py https://api-dot-isb-cgc.appspot.com/_ah/api/cohort_api/v1/patient_details?patient_barcode=TCGA-44-7670")
system(cmd, intern=T)

{
 "clinical_data": {
  "weiss_venous_invasion": "None",
  "primary_therapy_outcome_success": "Complete Remission/Response",
  "pregnancies": "None",
  "lymphatic_invasion": "None",
  "neoplasm_histologic_grade": "None",
  etc ...

```

For each barcode we get the list of associated files.

```
filePaths <- lapply(barcodes, function(bx) parseCurl(callCurl(bx)))
names(filePaths) <- barcodes

$`TCGA-NJ-A7XG-01A`
 [1] "gs://isb-cgc-open/file-path-currently-unavailable"
 [2] "gs://isb-cgc-open/tcga/luad/Genome_Wide_SNP_6/broad.mit.edu__snp_cnv/Level_3/NIZAM_p_TCGAb_406_07_08_09_10_N_GenomeWideSNP_6_B03_1486574.hg19.seg.txt"
 [3] "gs://isb-cgc-open/tcga/luad/Genome_Wide_SNP_6/broad.mit.edu__snp_cnv/Level_3/NIZAM_p_TCGAb_406_07_08_09_10_N_GenomeWideSNP_6_B03_1486574.nocnv_hg19.seg.txt"
 [4] "gs://isb-cgc-open/tcga/luad/HumanMethylation450/jhu-usc.edu__methylation/Level_3/jhu-usc.edu_LUAD.HumanMethylation450.18.lvl-3.TCGA-NJ-A7XG-01A-12D-A398-05.txt"
etc....

```

Now to filter down the files to a particular platform of interest

```
platformFilter <- function(filelist, filterby) {
  filelist[str_detect(filelist, filterby)]
}
```

Then we apply the filter to each set of file paths

```
rnaseqFiles <- lapply(filePaths, function(fp) platformFilter(fp, "rsem.genes.normalized_results"))

methylFiles <- lapply(filePaths, function(fp) platformFilter(fp, "HumanMethylation450"))
```

Looking at the rnaSeq files ... we have paths to google buckets.

```
$`TCGA-86-8076-01A`
[1] "gs://isb-cgc-open/tcga/luad/IlluminaHiSeq_RNASeqV2/unc.edu__RNASeqV2/Level_3/unc.edu.c9cf5f5f-30de-4062-99c3-dcb9cfd55870.1325628.rsem.genes.normalized_results"

$`TCGA-95-7948-01A`
[1] "gs://isb-cgc-open/tcga/luad/IlluminaHiSeq_RNASeqV2/unc.edu__RNASeqV2/Level_3/unc.edu.a5e6acbf-7027-4837-85d7-aacfe707900d.1371894.rsem.genes.normalized_results"

```

If a given barcode doesn't have files associated with the particular platform,
the list element will be empty.

Now, one could use [GCSfuse](https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=1&ved=0ahUKEwi7mMPY5prKAhVQ4WMKHcidAdIQFggdMAA&url=https%3A%2F%2Fgithub.com%2Fgooglecloudplatform%2Fgcsfuse&usg=AFQjCNH9eZaZ6IlcwYDL_7jW_jy4hM25wQ)
to mount a google bucket like "isb-cgc-open/tcga/luad/HumanMethylation450/jhu-usc.edu__methylation/Level_3/"
and work as usual. Alternatively, one can download associated files and work locally.
