#! /usr/bin/env python

import os,sys,re
from collections import defaultdict
import json
import shutil
import sqlite3
import pandas as pd
from Bio import SeqIO


def create_connection(db_file):
    """ create a database connection to the SQLite database
        specified by db_file
    :param db_file: database file
    :return: Connection object or None
    """
    conn = None
    try:
        conn = sqlite3.connect(db_file , check_same_thread=True)
        return conn
    except Error as e:
        print(e)

    return conn


directory = sys.argv[1]
sample2files = defaultdict(set)

json_filename = directory+'/'+'assembly'+'/'+'info.json'
if not os.path.exists(json_filename) :
    sys.exit(json_filename+' does not exist')


db_filename = directory+'/'+'refinedBins'+'/'+'REFINEDBINS.db'
if not os.path.exists(db_filename) :
    sys.exit(db_filename+' does not exist')


with open(json_filename) as f:
    data = json.load(f)
project = data['project']
sample = data['sample']
contig_filename = data['anvio_contig_filename']
protein_filename = data['assembly_protein_filename']
anvio_protein_filename = data['anvio_protein_filename']




genomeDir = 'genomes'
if os.path.exists(genomeDir) :
    sys.exit(genomeDir+' already exists, remove it first')
os.mkdir(genomeDir)

scaffold2coverage_filename = 'scaffold2coverage.info'
if os.path.exists(scaffold2coverage_filename) :
    sys.exit(genomeDir+' already exists, remove it first')


outputScaffold = open(scaffold2coverage_filename,'w')
outputScaffold.write('scaffold'+'\t'+'genome'+'\t'+'length_refineM'+'\t'+'coverage_refineM'+'\t'+'coverage_anvio'+'\n')



scaffold2bin = dict()
scaffold2refineM_cov = dict()
scaffold2anvio_cov = dict()
scaffold2len = dict()

conn = sqlite3.connect(db_filename , check_same_thread=True)

bin2name = dict()
sql_select_query = 'SELECT anvio_id,name FROM bins'
df_tmp = pd.read_sql_query(sql_select_query, con=conn)#.set_index('scaffold_id')
for index, row in df_tmp.iterrows():
    if row['anvio_id'] == 'Unbinned' :
        continue
    bin2name[ row['anvio_id'] ] = row['name']

sql_select_query = 'SELECT scaffold_id , anvio_updated_id , anvio_id , refineM_length , refineM_coverage , anvio_coverage FROM scaffolds'
df_tmp = pd.read_sql_query(sql_select_query, con=conn)#.set_index('scaffold_id')
for index, row in df_tmp.iterrows():
    if row['anvio_updated_id'] == 'Unbinned' :
        continue
    scaffold2bin[ row['scaffold_id'] ] = row['anvio_updated_id']
    scaffold2refineM_cov[ row['scaffold_id']  ] = row['refineM_coverage']
    scaffold2anvio_cov[ row['scaffold_id'] ] = row['anvio_coverage']
    scaffold2len[ row['scaffold_id'] ] = row['refineM_length']
conn.close()
print('close connection')

bin2seqList = defaultdict(list)
for record in SeqIO.parse(contig_filename,'fasta') :
    if record.id in scaffold2bin :
        bin2seqList[ bin2name[scaffold2bin[record.id]] ].append(record)

for name,seqList in bin2seqList.items() :
    genome_filename = genomeDir+'/'+name+'.fna'
    print(name+'\t'+str(len(seqList))+'\t'+genome_filename)
    SeqIO.write(seqList,genome_filename,'fasta')
    for record in seqList :
        outputScaffold.write(str(record.id)+'\t'+name+'\t'+str(scaffold2len[ record.id ])+'\t'+str(scaffold2refineM_cov[ record.id ])+'\t'+str(scaffold2anvio_cov[ record.id ])+'\n')
outputScaffold.close()



