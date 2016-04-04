FROM phusion/baseimage:0.9.18

CMD ["/sbin/my_init"]

# Add the PostgreSQL PGP key to verify their Debian packages.
# It should be the same key as https://www.postgresql.org/media/keys/ACCC4CF8.asc
RUN apt-key adv --keyserver hkp://p80.pool.sks-keyservers.net:80 --recv-keys B97B0AFCAA1A47F044F244A07FCC7D46ACCC4CF8

RUN echo "deb http://apt.postgresql.org/pub/repos/apt/ precise-pgdg main" > /etc/apt/sources.list.d/pgdg.list

RUN apt-get update && apt-get install -y \
    python-software-properties \
    software-properties-common \
    postgresql-9.4 \
    postgresql-client-9.4 \
    postgresql-contrib-9.4

RUN mkdir /source

ADD sql/projections.sql /source

WORKDIR /source

USER postgres

RUN    /etc/init.d/postgresql start &&\
    psql --command "CREATE USER projections_admin WITH SUPERUSER PASSWORD 'projections_password';" &&\
    createdb -O projections_admin projections_admin && createdb -O projections_admin projections_database && \
    psql projections_database -f projections.sql

# Adjust PostgreSQL configuration so that remote connections to the
# database are possible.
RUN echo "host all  all    0.0.0.0/0  md5" >> /etc/postgresql/9.4/main/pg_hba.conf

# And add ``listen_addresses`` to ``/etc/postgresql/9.3/main/postgresql.conf``
RUN echo "listen_addresses='*'" >> /etc/postgresql/9.4/main/postgresql.conf

# Expose the PostgreSQL port
EXPOSE 5432

# Add VOLUMEs to allow backup of config, logs and databases
VOLUME  ["/etc/postgresql", "/var/log/postgresql", "/var/lib/postgresql"]

# Set the default command to run when starting the container
CMD ["/usr/lib/postgresql/9.4/bin/postgres", "-D", "/var/lib/postgresql/9.4/main", "-c", "config_file=/etc/postgresql/9.4/main/postgresql.conf"]