#!/usr/bin/env python3


# This is temporary script created to init database in order to use tests
# Database projections_database must be created first.
# TODO remove this script

import getpass

import psycopg2

db_connection = psycopg2.connect("dbname=projections_database user={user_name}".format(user_name=getpass.getuser()))

cursor = db_connection.cursor()

cursor.execute("SELECT relname FROM pg_class WHERE relname = 'projections_table' ")
projections_table_exists = cursor.fetchone()

# If no projections table found, create required tables
if not bool(projections_table_exists):
    cursor.execute("""
    CREATE TABLE projections_table (
        projection_name varchar PRIMARY KEY UNIQUE,
        mount_path varchar UNIQUE,
        projector_pid int
    );""")
    db_connection.commit()

# Checking if tree_table exists in pg_class table
cursor.execute("SELECT relname FROM pg_class WHERE relname = %s;", ('tree_table',))
tree_table_exists = cursor.fetchone()

# If no tree table found, create required tables
if not bool(tree_table_exists):
    cursor.execute("CREATE TABLE tree_table ("
                   "node_id serial PRIMARY KEY, "
                   "projection_name varchar REFERENCES projections_table(projection_name) ON DELETE CASCADE,"
                   "parent_id integer, "
                   "name varchar, "
                   "uri varchar, "
                   "path varchar[][], "
                   "type varchar, "
                   "meta_links varchar[][]"
                   ");")
    db_connection.commit()

# Checking if metadata_table exists in pg_class table
cursor.execute("SELECT relname FROM pg_class WHERE relname = %s;", ('metadata_table',))
metadata_table_exists = cursor.fetchone()

if not metadata_table_exists:
    cursor.execute("CREATE TABLE metadata_table ("
                   "meta_id serial PRIMARY KEY,"
                   "node_id integer REFERENCES tree_table(node_id) ON DELETE CASCADE,"
                   "parent_node_id integer REFERENCES tree_table(node_id) ON DELETE CASCADE,"
                   "meta_contents jsonb);")
    db_connection.commit()

# Checking if projections_attributes_table exists in pg_class table
cursor.execute("SELECT relname FROM pg_class WHERE relname = %s;", ('projections_attributes_table',))
projections_attributes_table_exists = cursor.fetchone()

if not projections_attributes_table_exists:
    cursor.execute("CREATE TABLE projections_attributes_table ("
                   "attr_id serial PRIMARY KEY,"
                   "node_id integer REFERENCES tree_table(node_id) ON DELETE CASCADE, "
                   "st_atime int, "
                   "st_mtime int, "
                   "st_ctime int, "
                   "st_size int, "
                   "st_mode varchar, "
                   "st_nlink int, "
                   "st_ino int);")
    db_connection.commit()
