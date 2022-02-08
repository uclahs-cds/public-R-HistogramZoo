FROM r-base
RUN apt-get update
RUN apt-get -y install libcurl4-openssl-dev libssl-dev libxml2 libxml2-dev libicu-dev \
                       libpq-dev libmysqlclient-dev curl \
                       libcairo2-dev libxt-dev libproj-dev pandoc
RUN Rscript -e 'install.packages("devtools")'
RUN Rscript -e 'install.packages("BiocManager")'
RUN Rscript -e 'BiocManager::install(c("GenomicRanges", "GenomicFeatures", "metap", "S4Vectors", "AnnotationDbi"))'
RUN chmod 777 /usr/local/lib/R/site-library
COPY . /usr/local/src/ConsensusPeaks
WORKDIR /usr/local/src/ConsensusPeaks
RUN Rscript -e 'devtools::install_local(dependencies = T, force = T)'
