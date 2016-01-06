# Accessing low level *controlled* data with R/Bioconductor

To access controlled data, we need membership to the google group called isb-cgc-cntl (compared to isb-cgc-open).
To join the group, we need to verify our dbGaP credentials:

(a) sign in to isb-cgc.appspot.com
(b) After signing in, click on the circle/picture next to your name in the upper-right corner of the app
(c) in the middle, in the "Data Access" page, click on the "Associate with eRA Commons Account" line
(d) this should redirect you to the NIH secure iTrust site
(e) enter your username and password
(f) you will magically come back to the ISB-CGC app page with a blue warning box about accessing TCGA controlled data

After clicking on the circle/picture thing again,
you should see a message that tells you how much time you're authorized for
(24 hours from when you authenticated thru NIH)

Now your google identity should be on the isb-cgc-cntl google group which has access to the controlled-access CEL files so try, for example:

gs://62f2c827-mock-mock-mock-1cde698a4f77/tcga/acc/Genome_Wide_SNP_6/broad.mit.edu__snp_cnv/Level_1/AQUAE_p_TCGA_112_304_b2_N_GenomeWideSNP_6_F04_1348258.CEL

(which is not *really* controlled TCGA data but just "mock")

and:

gs://62f2c827-93cc-4ca7-a90f-1cde698a4f77/tcga/acc/Genome_Wide_SNP_6/broad.mit.edu__snp_cnv/Level_1/AQUAE_p_TCGA_112_304_b2_N_GenomeWideSNP_6_A01_1348356.CEL


## Starting Docker and mounting data

So let's access some data!

```
Starting machine default...
(default) OUT | Starting VM...
...
                        ##         .
                  ## ## ##        ==
               ## ## ## ## ##    ===
           /"""""""""""""""""\___/ ===
      ~~~ {~~ ~~~~ ~~~ ~~~~ ~~~ ~ /  ===- ~~~
           \______ o           __/
             \    \         __/
              \____\_______/


docker is configured to use the default machine with IP xxx.yyy.zz.www
For help getting started, check out the docs at https://docs.docker.com
```

Now we can start up the docker image, assuming you have already pulled the
bioconductor microarray container (docker pull bioconductor/release_microarray).

```
docker run -ti --privileged bioconductor/release_microarray /bin/bash
curl https://sdk.cloud.google.com | bash
bash
gcloud init
```


OK, now, we need to be able to access data that sitting in a controlled access
bucket. One way is to use gsutil to copy the data to our local drive.


```
gsutil cp gs://62f2c827-93cc-4ca7-a90f-1cde698a4f77/tcga/acc/\
Genome_Wide_SNP_6/broad.mit.edu__snp_cnv/Level_1/\
AQUAE_p_TCGA_112_304_b2_N_GenomeWideSNP_6_A01_1348356.CEL .
```

We can also use GCSFuse to mount the bucket, and work with it interactively.
Different from working with the open access bucket, we need to use the
"--implicit-dirs" flag to get proper directory listings.


```
apt-get install fuse curl daemon
curl -L -O https://github.com/GoogleCloudPlatform/gcsfuse/releases/download/v0.12.0/gcsfuse_0.12.0_amd64.deb
sudo dpkg --install gcsfuse_0.12.0_amd64.deb
mkdir /media/dat
daemon -- gcsfuse --implicit-dirs 62f2c827-93cc-4ca7-a90f-1cde698a4f77 /media/dat
```


Great! We got it! Now we can use Bioconductor workflows to process the data.


## Locating data using the GCG web service

First we need to find the files. Let's suppose we are interested in a particular
sample ID (Or set of IDs). We will use the API endpoint "datafilenamekey_list" to
get file paths given a TCGA barcode.

But, as Tyrell says in Bladerunner, "I want to see a negative before I provide you with a positive."
So let's do that.

```
R
install.packages("httr")
library(httr)
r <- GET("https://mvm-dot-isb-cgc.appspot.com/_ah/api/cohort_api/v1/datafilenamekey_list?sample_barcode=TCGA-Bad-Data")
r
```

And the response is ...

```
Response [https://mvm-dot-isb-cgc.appspot.com/_ah/api/cohort_api/v1/datafilenamekey_list?sample_barcode=TCGA-Bad-ID]
  Date: 2015-12-16 22:51
  Status: 200
  Content-Type: application/json; charset=UTF-8
  Size: 125 B
{
 "count": "0",
 "kind": "cohort_api#cohortsItem",
 "etag": "\"YSI636SyFAiMSbxaGP932XeRcUk/3_-4ZVfBVVlUTjAMHiYmkJAwlE0\""
}
```

The "count" is the number of files associated with the Barcode, and here, with a bad barcode,
we get zero files asscociated. So what about a valid id?

```
r <- GET("https://isb-cgc.appspot.com/_ah/api/cohort_api/v1/datafilenamekey_list?sample_barcode=TCGA-W5-AA31-01A")
r
```
Here's the response:

```
Response [https://isb-cgc.appspot.com/_ah/api/cohort_api/v1/datafilenamekey_list?sample_barcode=TCGA-W5-AA31-01A]
  Date: 2015-12-22 07:04
  Status: 200
  Content-Type: application/json; charset=UTF-8
  Size: 2.1 kB
{
 "count": "13",
 "datafilenamekeys": [
  "gs://isb-cgc-open/file-path-currently-unavailable",
  "gs://isb-cgc-open/tcga/chol/Genome_Wide_SNP_6/broad.mit.edu__snp_cnv/Level...
  "gs://isb-cgc-open/tcga/chol/Genome_Wide_SNP_6/broad.mit.edu__snp_cnv/Level...
  "gs://isb-cgc-open/tcga/chol/HumanMethylation450/jhu-usc.edu__methylation/L...
  "gs://isb-cgc-open/tcga/chol/IlluminaHiSeq_miRNASeq/bcgsc.ca__miRNASeq/Leve...
  "gs://isb-cgc-open/tcga/chol/IlluminaHiSeq_miRNASeq/bcgsc.ca__miRNASeq/Leve...
  "gs://isb-cgc-open/tcga/chol/IlluminaHiSeq_RNASeqV2/unc.edu__RNASeqV2/Level...
...
```

Great! We have 20 files associated with this barcode. Among the data we have
SNP 6.0, methylation array data, VCFs (DNAseq), RNAseq, miRNAseq, and RPPA data. Let's
try accessing one of the CEL files.

## Processing data with Bioconductor

To work with the Affy SNP 6.0 CEL files, we're going to use the Bioconductor
package "oligo". To properly read the files, we're going to have to download
a large annotation package.

```
source("https://bioconductor.org/biocLite.R")
biocLite("pd.genomewidesnp.6")
library(oligo)
library(pd.genomewidesnp.6)

celfiles <- list.files("/media/dat/ccle/SNP_Arrays/", pattern=".CEL")
celpaths <- sapply(celfiles[1:3], function(a) paste("/media/cntl/tcga/acc/Genome_Wide_SNP_6/broad.mit.edu__snp_cnv/Level_1/", a, sep=""))

rawData <- read.celfiles(celpaths)

pdf("image_cel1.pdf")
image(rawData, which=1, transfo=log2)
dev.off()
```

Success!!  Now we need to get that image out of the docker. After we exit from
R, we can see the file. Without logging out of the docker container, I'm going
to open another terminal window, and copy out that file.

```
eval "$(docker-machine env default)"
docker ps # to get the container ID
docker cp 1cf74ce69172:image_cel1.pdf .
```

So, we mounted a bucket in a docker container. Read a raw PROTECTED CEL file from that
bucket, and using a bioconductor package, produced an image of that microarray.
And finally copied it out to our local system.
