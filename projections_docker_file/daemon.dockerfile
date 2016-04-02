FROM phusion/baseimage:0.9.18

CMD ["/sbin/my_init"]

ENV LANG en_US.UTF-8

RUN apt-get update && apt-get install -y \
    gcc \
    python3 python3-dev python3-pip \
    python-software-properties \
    wget \
    git \
    tar \
    libpq-dev \
    libfuse-dev

COPY projections-core.tar.gz /

RUN tar -xzf projections-core.tar.gz

WORKDIR /projections-core

RUN pip3 install -r requirements.txt

RUN ["python3", "/projections-core/projections_daemon.py", "-start"]
