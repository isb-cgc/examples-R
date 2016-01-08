# Accessing low level *controlled* data

To access controlled data, you need to be a member of the isb-cgc-cntl google group (compared to isb-cgc-open).
To join the group, you need to verify your dbGaP credentials:

1. Sign in to isb-cgc.appspot.com
2. After signing in, click on the circle/picture next to your name in the upper-right corner of the app
3. In the middle, in the "Data Access" page, click on the "Associate with eRA Commons Account" line
4. This should redirect you to the NIH secure iTrust site
5. Enter your username and password
6. You will magically come back to the ISB-CGC app page with a blue warning box about accessing TCGA controlled data

After clicking on the circle/picture thing again, you should see a message that tells you how much time you're authorized for (24 hours from when you authenticated thru NIH)

## Starting Docker and getting set up

Let's get started! You might want to take a look at the Processing_Raw_Data_With_Bioconductor.md
file for more information about working with Docker.

The first thing we need to do is start up docker.

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

An easy way to move data in and out of docker, involves mounting a data volume.
Check the documentation for details [docker data volumes](https://docs.docker.com/engine/userguide/dockervolumes/).

```
#<< run on the command line >>#
mkdir /my_docker_workspace
docker run -ti -v=/my_docker_workspace:workspace --privileged bioconductor/release_microarray /bin/bash
```

Then, inside docker, the directory /workspace will mirror what you have in
your userspace director. This makes it easy to move files in and out.
Now, we'll get authorized.
Notice we're downloading the 'raw' files from github. We're going to have to
authorize via gcloud, and then isb_auth.py.

```
#<< run each of these on the command line >>#

# downloads the gcloud application
curl https://sdk.cloud.google.com | bash

# need to start a new bash session to "see" gcloud
bash

# authorize the session
gcloud init

# the isb authorization scrips need the oauth2client library
pip install --upgrade oauth2client

# download and run the isb authorization script.
wget https://raw.githubusercontent.com/isb-cgc/ISB-CGC-Webapp/master/scripts/isb_auth.py
python isb_auth.py --noauth_local_webserver

# then download the isb script to access file lists.
wget https://raw.githubusercontent.com/isb-cgc/ISB-CGC-Webapp/master/scripts/isb_curl.py
```

The gcloud init and the isb_auth.py script will both give you
long web links to feed into a browser, which returns verification codes from google.
At the moment, we need to verify once to get restricted file lists from ISB, and
once to get to the restricted google buckets containing the data. In the future,
this will likely be a single step.

## Locating data using the GCG web service

Let's open up R!
For an example, we're going to look up the files associated with a given barcode,
TCGA-W5-AA31-01A.

```
library(stringr)
x <- system("python isb_curl.py https://isb-cgc.appspot.com/_ah/api/cohort_api/v1/datafilenamekey_list?sample_barcode=TCGA-W5-AA31-01A", intern=T)
```

We're going to need a little function to parse through what returns.

```
parseCurl <- function(x) {
  y <- str_trim(strsplit(x, ",")[[1]])
  idx <- str_detect(pattern="^u'gs://", string=y)  # lines containing files
  z <- str_split(y[idx], "\'")
  unlist(lapply(z, function(zi) zi[2]))
}

fileList <- parseCurl(x)
cmd <- paste0("gsutil cp ", fileList[3], " /media/dat")
system(cmd)
```

Great! We have 22 files associated with this barcode. Among the data we have
SNP 6.0, methylation array data, VCFs (DNAseq), RNAseq, miRNAseq, and RPPA data.
We downloaded one of the CEL files to /media/dat... now we can process it!

## Processing data with Bioconductor

To work with the Affy SNP 6.0 CEL files, we're going to use the Bioconductor
package "oligo". To properly read the files, we're going to have to download
a large annotation package.

```
source("https://bioconductor.org/biocLite.R")
biocLite("pd.genomewidesnp.6")
library(oligo)
library(pd.genomewidesnp.6)

celfiles <- list.files("/media/dat", pattern=".CEL")

rawData <- read.celfiles(celfiles)

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

So, gathered all the files, protected and open, associated with a barcode.
Downloaded a raw PROTECTED CEL file, and using a bioconductor package,
produced an image of that microarray. And finally copied it out to our local system.
