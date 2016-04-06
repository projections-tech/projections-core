# This is sandbox project for FUSE based projections and metadata system

##Requirements
Any modern Linux machine with steady internet connection will suffice.

##Installation

In order to run tests create database named 'projections_database' using following command:

'''
$ createdb projections_database
'''

Then run database setup SQL script:
'''
$ psql projections_database -f sql/projections.sql
'''

In order to run SRA, system specific SRA-toolkit (for example sratoolkit.2.5.2-ubuntu64) must be placed in project root folder. SRA-toolkit can be found [here](http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software)
