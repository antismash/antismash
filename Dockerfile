# antiSMASH 5 container with a snapshot of the development tree
# VERSION 0.0.1
FROM antismash/base5-nonfree:latest
LABEL maintainer="Kai Blin <kblin@biosustain.dtu.dk>"

# Python and Docker are not getting along encoding-wise
ENV LANG C.UTF-8

# Install git
RUN apt-get update && apt-get install -y git && apt-get clean -y && apt-get autoclean -y && apt-get autoremove -y && rm -rf /var/lib/apt/lists/*

# Grab antiSMASH
COPY . /antismash

ADD docker/instance.cfg /antismash/antismash/config

RUN HARDCODE_ANTISMASH_GIT_VERSION=1 pip3 install /antismash && python3 -c "import antismash; antismash.config.build_config(['--databases', 'mounted_at_runtime'], modules=antismash.get_all_modules()); antismash.main.prepare_module_data()"

RUN mkdir /matplotlib && MPLCONFIGDIR=/matplotlib python3 -c "import matplotlib.pyplot as plt" && chmod -R a+rw /matplotlib

ADD docker/run /usr/local/bin/run

VOLUME ["/input", "/output", "/databases"]
WORKDIR /output

ENTRYPOINT ["/usr/local/bin/run"]
