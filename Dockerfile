# Base image https://hub.docker.com/u/rocker/
FROM rocker/r-base:latest

## install debian packages
RUN apt-get update -qq && apt-get -y --no-install-recommends install \
libxml2-dev \
libcairo2-dev \
libsqlite3-dev \
libmariadbd-dev \
libpq-dev \
libssh2-1-dev \
unixodbc-dev \
libcurl4-openssl-dev \
libssl-dev

## create directories (Mirror directory structure)
RUN mkdir -p /resources
RUN mkdir -p /src
RUN mkdir -p /output

## copy files
COPY /src/install_packages.R /src/install_packages.R
COPY /src/sambamba_exon_coverage.R /src/sambamba_exon_coverage.R 

## install R-packages
RUN Rscript /src/install_packages.R

