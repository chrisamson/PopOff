# build the Docker image from the base image 'openanalytics/r-base'
# this is an Ubuntu 16.04 LTS with a recent R version.
# this image is available on Docker hub at https://hub.docker.com/r/openanalytics/r-base/
FROM openanalytics/r-ver:4.1.3

LABEL maintainer "Christopher Samson <chrisamson10@gmail.com>"

# system libraries of general use
RUN apt-get update && apt-get install --no-install-recommends -y \
    pandoc \
    pandoc-citeproc \
    libcurl4-gnutls-dev \
    libcairo2-dev \
    libxt-dev \
    libssl-dev \
    libssh2-1-dev \
    libssl1.1 \
    cmake \
    && rm -rf /var/lib/apt/lists/*

# install basic shiny functionality to R
RUN R -e "install.packages(c('RSQLite','shiny','DT','purr','ggplot2','scales','ggpubr','tidyr','dplyr','plyr','grid','gridExtra'))" #, repos='https://cloud.r-project.org/')"

COPY POPOFF.R /root/POPOFF.R
COPY Rprofile.site /usr/local/lib/R/etc/
COPY POPOFF_DB /root/POPOFF_DB
WORKDIR /root

EXPOSE 3838

CMD ["Rscript","/root/POPOFF.R"]
#CMD ["R", "-q", "-e", "shiny::runApp('/root/POPOFF')"]