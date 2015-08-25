Docker Image Scripts
=============================================

[Bioconductor](http://www.bioconductor.org/) provides an excellent set of [docker images](http://www.bioconductor.org/help/docker/) containing R, RStudioServer, and the sets of Bioconductor packages appropriate for certain use cases.

Here we have a version of the Bioconductor `bioconductor/devel_base` which adds on R packages containing training material such as [ISB CGC Examples](https://github.com/isb-cgc/examples-R), [gcloud](https://cloud.google.com/sdk/gcloud/), and some configuration to .Rprofile.

Deploying a new version to the public image repository
------------------------------------------------------
Get the latest Dockerfile, etc., from this repository via a `git clone` or `git pull` and then:

(0) Make sure your build machine has the latest Bioconductor Docker image.
```
sudo docker pull bioconductor/devel_base
```

(1) Build the image.
```
sudo docker build -t gcr.io/isb-cgc/r-examples:0.01 .
```

(2) Push the new version to the public image repository.  *Always specify a tag.*
```
sudo gcloud docker push gcr.io/isb-cgc/r-examples:0.01
```

(3) Also tag the new version as 'latest'.  *Always explicity mark as 'latest' a particular tagged version.*
```
sudo docker tag  gcr.io/isb-cgc/r-examples:0.01 gcr.io/isb-cgc/r-examples:latest
```

(4) And push 'latest'. (This will be really quick since its just updating metadata about 'latest'.)
```
sudo gcloud docker push gcr.io/isb-cgc/r-examples:0.01:latest 
```
