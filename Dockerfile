FROM bioconductor/bioconductor_docker:latest
RUN apt-get update && apt-get install libcurl4-openssl-dev libxml2-dev
RUN wget 'http://downloads.sourceforge.net/project/sbml/libsbml/5.11.4/stable/libSBML-5.11.4-core-src.tar.gz?r=&ts=1435978044' -O /tmp/libsbml.tar.gz && \
 tar -xf /tmp/libsbml.tar.gz && \
 cd libsbml-5.11.4 && ./configure --prefix=/usr/local >/dev/null && sudo make install >/dev/null && cd .. && sudo rm -rf libsbml-5.11.4

RUN R -e "BiocManager::install('NetPathMiner', ask=F)"
