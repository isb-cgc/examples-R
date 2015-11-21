Docker Image Scripts
=============================================

[Bioconductor](http://www.bioconductor.org/) provides an excellent set of [docker images](http://www.bioconductor.org/help/docker/) containing R, RStudioServer, and the sets of Bioconductor packages appropriate for certain use cases.

Here we have a version of the Bioconductor `bioconductor/release_core` which adds on R packages containing training material such as [ISB CGC Examples](https://github.com/isb-cgc/examples-R), [gcloud](https://cloud.google.com/sdk/gcloud/), and some configuration to .Rprofile.

Deploying a new version to the public image repository
------------------------------------------------------

The following instructions assume that Docker is installed and the current user can run docker commands.  For help with that one-time initial setup, see https://cloud.google.com/container-registry/docs/#install_docker.

(1) Get the latest Dockerfile, etc., from this repository via a `git clone https://github.com/isb-cgc/examples-R.git` or `git pull`.

(2) Change into the directory where the Dockerfile reside.
```
cd examples-R/inst/docker
```

(3) Make sure your build machine has the latest Bioconductor Docker image.
```
docker pull bioconductor/release_core
```

(4) Build the image using a tag indicating today's date. *Always specify a tag.*
```
docker build -t b.gcr.io/isb-cgc-public-docker-images/r-examples:2015-10-30 .
```

(5) Push the new version to the public image repository.
```
gcloud docker push b.gcr.io/isb-cgc-public-docker-images/r-examples:2015-10-30
```

(6) Also tag the new version as 'latest'.  *Always explicity mark as 'latest' a particular tagged version.*
```
docker tag -f b.gcr.io/isb-cgc-public-docker-images/r-examples:2015-10-30 \
  b.gcr.io/isb-cgc-public-docker-images/r-examples:latest
```

(7) And push 'latest'. (This will be really quick since its just updating metadata about 'latest'.)
```
gcloud docker push b.gcr.io/isb-cgc-public-docker-images/r-examples:latest 
```
