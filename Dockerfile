FROM ubuntu:19.10 as builder

USER root

RUN apt-get -yq update
RUN apt-get install -yq --no-install-recommends \
locales \
g++ \
gcc \
r-base \
r-cran-devtools \
make \
libcurl4-openssl-dev \
libgit2-dev \
zlib1g-dev \
libxml2-dev \
libssl-dev \
libbz2-dev \
liblzma-dev \
gfortran


RUN locale-gen en_US.UTF-8
RUN update-locale LANG=en_US.UTF-8

ENV OPT /opt/local
ENV PATH $OPT/bin:$PATH
ENV R_LIBS $OPT/R-lib
ENV R_LIBS_USER $R_LIBS
ENV LD_LIBRARY_PATH $OPT/lib
ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8

RUN mkdir -p $R_LIBS
RUN mkdir -p $LD_LIBRARY_PATH

# build tools from other repos
ADD R/libInstall.R build/
ADD . deepSNV/
RUN Rscript build/libInstall.R


FROM ubuntu:19.10

RUN apt-get -yq update
RUN apt-get install -yq --no-install-recommends \
locales \
curl \
r-base \
r-cran-devtools \
zlib1g \
libxml2 \
bzip2 \
unattended-upgrades && \
unattended-upgrade -d -v && \
apt-get remove -yq unattended-upgrades && \
apt-get autoremove -yq

RUN locale-gen en_US.UTF-8
RUN update-locale LANG=en_US.UTF-8

ENV OPT /opt/local
ENV PATH $OPT/bin:$PATH
ENV R_LIBS $OPT/R-lib
ENV R_LIBS_USER $R_LIBS
ENV LD_LIBRARY_PATH $OPT/lib
ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8

RUN mkdir -p $OPT
COPY --from=builder $OPT $OPT

WORKDIR /home/ubuntu

CMD ["/bin/bash"]
