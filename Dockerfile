# antiSMASH container with a snapshot of the development tree
# VERSION 0.0.2
FROM docker.io/antismash/base:latest
LABEL maintainer="Kai Blin <kblin@biosustain.dtu.dk>"

# Python and Docker are not getting along encoding-wise
ENV LANG C.UTF-8

# set up antiSMASH deb repo
ADD https://dl.secondarymetabolites.org/antismash-stretch.list /etc/apt/sources.list.d/antismash.list
ADD https://dl.secondarymetabolites.org/antismash.asc /tmp/
RUN apt-key add /tmp/antismash.asc

# Install git and meme-suite
RUN apt-get update && apt-get install -y git meme-suite && apt-get clean -y && apt-get autoclean -y && apt-get autoremove -y && rm -rf /var/lib/apt/lists/*

# Grab antiSMASH
COPY . /antismash

ADD docker/instance.cfg /antismash/antismash/config

RUN HARDCODE_ANTISMASH_GIT_VERSION=1 pip3 install /antismash --break-system-packages && python3 -c "import antismash; antismash.config.build_config(['--databases', 'mounted_at_runtime'], modules=antismash.get_all_modules()); antismash.main.prepare_module_data()"

RUN mkdir /matplotlib && MPLCONFIGDIR=/matplotlib python3 -c "import matplotlib.pyplot as plt" && chmod -R a+rw /matplotlib

ADD docker/run /usr/local/bin/run

VOLUME ["/input", "/output", "/databases"]
WORKDIR /output

ENTRYPOINT ["/usr/local/bin/run"]
