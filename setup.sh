#!/usr/bin/env bash

# This script performs projection database setup.
if psql -lqt | cut -d \| -f 1 | grep -qw projections_database; then
    echo 'Database exists!'
else
    # Create database and initialize it
    createdb projections_database && psql projections_database -f projection_database_inititalization.sql
fi
