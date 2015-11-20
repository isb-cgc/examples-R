### Start from Bioconductor release core ######
FROM bioconductor/release_core
MAINTAINER deflaux@google.com

### Add training material packages #####################################
# Modelled after https://github.com/Bioconductor/bioc_docker/blob/master/out/devel_sequencing/Dockerfile
ADD .Rprofile /home/rstudio/

ADD installPackages.R /tmp/
RUN cd /tmp && \
    R -f /tmp/installPackages.R

### Install gcloud ##############################################
# Modelled after https://github.com/GoogleCloudPlatform/cloud-sdk-docker/blob/master/Dockerfile
RUN apt-get update && \
    apt-get install -y -qq --no-install-recommends \
      openjdk-7-jre-headless \
      openssh-client \
      php5-cgi \
      php5-cli \
      php5-mysql \
      python \
      python-openssl \
      unzip \
      wget && \
    apt-get clean
RUN wget https://dl.google.com/dl/cloudsdk/release/google-cloud-sdk.zip && \
  unzip google-cloud-sdk.zip && \
  rm google-cloud-sdk.zip
ENV CLOUDSDK_PYTHON_SITEPACKAGES 1
RUN google-cloud-sdk/install.sh --usage-reporting=true \
  --path-update=true \
  --bash-completion=true \
  --rc-path=/.bashrc \
  --disable-installation-options
RUN google-cloud-sdk/bin/gcloud --quiet components update \
  pkg-go \
  pkg-python \
  pkg-java \
  preview \
  alpha \
  beta \
  app
RUN google-cloud-sdk/bin/gcloud --quiet config set component_manager/disable_update_check true
RUN mkdir /.ssh
ENV PATH /google-cloud-sdk/bin:$PATH
