#! /usr/bin/env python


import os,sys,re
import sqlite3
from sqlite3 import Error


def create_connection(db_file):
    """ create a database connection to a SQLite database """
    conn = None
    try:
        conn = sqlite3.connect(db_file)
        print(sqlite3.version)
    except Error as e:
        print(e)
    finally:
        if conn:
            conn.close()


if __name__ == '__main__':
    directory = '/env/cns/proj/projet_CSD/scratch/assemblies/Esil_F_AA1'
    db_filename = directory+'/'+'refinedBins'+'/'+'REFINEDBINS.db'
    create_connection(db_filename)
