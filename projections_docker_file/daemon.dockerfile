FROM phusion/baseimage:0.9.18

CMD ["/sbin/my_init"]

ENV LANG en_US.UTF-8

# Add the PostgreSQL PGP key to verify their Debian packages.
# It should be the same key as https://www.postgresql.org/media/keys/ACCC4CF8.asc
RUN apt-key adv --keyserver hkp://p80.pool.sks-keyservers.net:80 --recv-keys B97B0AFCAA1A47F044F244A07FCC7D46ACCC4CF8

RUN echo "deb http://apt.postgresql.org/pub/repos/apt/ precise-pgdg main" > /etc/apt/sources.list.d/pgdg.list

RUN apt-get update && apt-get install -y \
    gcc \
    python3 python3-dev python3-pip \
    python-software-properties \
    wget \
    git \
    tar \
    libpq-dev \
    libfuse-dev

RUN mkdir /source

COPY projections-core.tar.gz /source

WORKDIR /source

RUN tar -xzf projections-core.tar.gz

RUN pip3 install -r projections-core/requirements.txt

RUN /bin/bash

