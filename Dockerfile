# DO NOT CHANGE
from 812206152185.dkr.ecr.us-west-2.amazonaws.com/latch-base:fe0b-main

workdir /tmp/docker-build/work/

shell [ \
    "/usr/bin/env", "bash", \
    "-o", "errexit", \
    "-o", "pipefail", \
    "-o", "nounset", \
    "-o", "verbose", \
    "-o", "errtrace", \
    "-O", "inherit_errexit", \
    "-O", "shift_verbose", \
    "-c" \
]
env TZ='Etc/UTC'
env LANG='en_US.UTF-8'

arg DEBIAN_FRONTEND=noninteractive

# Latch SDK
# DO NOT REMOVE
run pip install latch==2.50.3
run mkdir /opt/latch

RUN apt-get update -y && \
apt-get install -y autoconf curl zip unzip wget gcc git make libbz2-dev zlib1g-dev pip \
libncurses5-dev libncursesw5-dev liblzma-dev libtool autoconf build-essential pkg-config automake tcsh

RUN pip install pandas
# Copy workflow data (use .dockerignore to skip files)
RUN wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh && \
    bash Mambaforge-Linux-x86_64.sh -b -p /mambaforge && \
    rm Mambaforge-Linux-x86_64.sh
ENV PATH="/mambaforge/bin:$PATH"
RUN mamba create -n crispresso2_env -c conda-forge -c bioconda crispresso2

copy . /root/


# Latch workflow registration metadata
# DO NOT CHANGE
arg tag
# DO NOT CHANGE
env FLYTE_INTERNAL_IMAGE $tag

workdir /root
