# Base image https://hub.docker.com/u/rocker/
FROM rocker/verse:latest

## create directories (Mirror directory structure, /src already exists)
RUN mkdir -p /resources
RUN mkdir -p /output

## copy files
COPY /src/install_packages.R /src/install_packages.R
COPY /src/sambamba_exon_coverage.R /src/sambamba_exon_coverage.R 

## install R-packages
RUN Rscript /src/install_packages.R

