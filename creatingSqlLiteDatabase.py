#! /usr/bin/env python


import os,sys,re
import sqlite3
from sqlite3 import Error
from collections import defaultdict


def create_connection(db_file):
    """ create a database connection to the SQLite database
        specified by db_file
    :param db_file: database file
    :return: Connection object or None
    """
    conn = None
    try:
        conn = sqlite3.connect(db_file)
        return conn
    except Error as e:
        print(e)

    return conn


def create_table(conn, create_table_sql):
    """ create a table from the create_table_sql statement
    :param conn: Connection object
    :param create_table_sql: a CREATE TABLE statement
    :return:
    """
    print(create_table_sql)
    try:
        c = conn.cursor()
        c.execute(create_table_sql)
    except Error as e:
        print(e)



def suggestedNameFunction(anvio_lineage,gtdb_lineage,project,sample,binName2count) :
    if gtdb_lineage != 'NULL' :
        liste = gtdb_lineage.split(';')
        if re.match(r's\_\_',liste[-1]) :
            if liste[-1] == 's__' : # if gtdb species name is empty, look for genus, family, order, class, phylum, domain name...
                name = 'Unknown'
                #print(reversed(liste))
                for taxa in reversed(liste) :
                    if taxa.split('__')[1] != '' :
                        name = taxa.split('__')[1]
                        break
                    else:
                        continue
            else: # if has a gtdb species name then concatenate genus and species name
                name = liste[-1].split('__')[1].replace(' ','_')
        else: # no species s__ field...
            name = liste[-1].split('__')[1]
        name = name.replace('-','_')
    else:
        name = anvio_lineage.replace(' ','_')
        name = name.replace('-','_')
        
    binName = name+'__'+project+'__'+sample
    # check name redundancy
    if binName in binName2count :
        nb = str(binName2count[binName])
        binName2count[binName] += 1
        binName = name+'_'+nb+'__'+project+'__'+sample
    else:
        binName2count[binName] = 1
    return binName


def create_bin(conn,binList) :
    """
    Create a new bin into the bins table
    :param conn:
    :param bin:
    :return: bin id
    """
    #print(binList)
    #print(len(binList))
    sql = ''' INSERT INTO bins (anvio_id,name,anvio_coverage,anvio_length,anvio_gc,anvio_contig_nb,anvio_N50,anvio_completeness,anvio_contamination,anvio_taxonomy,gtdb_taxonomy,gtdb_fastani_reference,gtdb_classification_method,checkM_nb_predicted_genes,checkM_translation_table,checkM_coding_density,checkM_completeness,checkM_contamination,anvio_bin)
              VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?) '''
    cur = conn.cursor()
    cur.execute(sql, binList)
    conn.commit()
    #print(cur.lastrowid)
    return cur.lastrowid

def populate_bins_table(conn, bins_summary_filename, sample,project) :
    binName2count = defaultdict(int)
    file = open(bins_summary_filename,'r')
    headerList = next(file).rstrip().split('\t')
    for line in file :
        line  = line.rstrip()
        liste = line.split('\t')
        anvio_id = liste[0]
        anvio_taxonomy = liste[1]
        anvio_coverage = float(liste[2])
        anvio_length = int(liste[3])
        anvio_contig_nb = int(liste[4])
        anvio_N50 = int(liste[5])
        anvio_gc = float(liste[6])/100.0
        anvio_completeness = float(liste[7])/100.0
        anvio_contamination = float(liste[8])/100.0

        if liste[9] == 'Na' :
            checkM_nb_predicted_genes = 'NULL'
        else:
            checkM_nb_predicted_genes = int(liste[9])

        if liste[10] == 'Na' :
            checkM_translation_table = 'NULL'
        else:
            checkM_translation_table = int(liste[10])

        if liste[11] == 'Na' :
            checkM_coding_density = 'NULL'
        else:
            checkM_coding_density = float(liste[11])

        if liste[12] == 'Na' :
            checkM_completeness = 'NULL'
        else:
            checkM_completeness = float(liste[12])/100.0

        if liste[13] == 'Na' :
            checkM_contamination = 'NULL'
        else:
            checkM_contamination = float(liste[13])/100.0

        if liste[14] == 'Na' :
            gtdb_taxonomy = 'NULL'
        else:
            gtdb_taxonomy = liste[14]

        if liste[15] == 'Na' :
            gtdb_fastani_reference = 'NULL'
        else:
            gtdb_fastani_reference = liste[15]

        if liste[16] == 'Na' :
            gtdb_classification_method = 'NULL'
        else:
            gtdb_classification_method = liste[16]


        name = suggestedNameFunction(anvio_taxonomy,gtdb_taxonomy,project,sample,binName2count)

        row_id = create_bin(conn,[anvio_id,name,anvio_coverage,anvio_length,anvio_gc,anvio_contig_nb,anvio_N50,anvio_completeness,anvio_contamination,anvio_taxonomy,gtdb_taxonomy,gtdb_fastani_reference,gtdb_classification_method,checkM_nb_predicted_genes,checkM_translation_table,checkM_coding_density,checkM_completeness,checkM_contamination , 1])
        #print(row_id)
    file.close()


def create_scaffold(conn,scaffoldList) :
    """
    Create a new scaffold into the scaffolds table
    :param conn:
    :param bin:
    :return: bin id
    """
    #print(scaffoldList)
    #print(len(scaffoldList))
    sql = ''' INSERT INTO scaffolds (scaffold_id , anvio_id , anvio_updated_id , anvio_gc , anvio_nb_splits , anvio_coverage , anvio_taxonomy , anvio_length , refineM_outlier , refineM_length , refineM_coverage , refineM_gc , refineM_nb_genes  ,refineM_coding_length , refineM_domain , refineM_phylum , refineM_class , refineM_order , refineM_family ,refineM_genus , refineM_species )
              VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?) '''
    cur = conn.cursor()
    cur.execute(sql, scaffoldList)
    conn.commit()
    #print(cur.lastrowid)
    return cur.lastrowid


def populate_scaffolds_table(conn, directory) :
    binFilenameSet = set()
    filenameSet = set( [ 'Collection.xlsx' , 'Anvio_summary.tsv' , 'CheckM.tsv' , 'GTDBtk.tsv' , 'Bins_summary.tsv' , 'Sample_summary.tsv' , 'partial_scaffolds.tsv' ] )
    for root, dirs, files in os.walk(directory):
        for filename in files :
            if filename in filenameSet :
                continue
            else:
                binFilenameSet.add(root+'/'+filename)


    data_bin = list()
    for filename in binFilenameSet :
        lengthSet = set()
        anvio_id = os.path.basename(filename.replace('.tsv',''))
        file = open(filename,'r')
        headerList = next(file).rstrip().split('\t')
        for line in file :
            line = line.rstrip()
            liste = line.split('\t')
        
            scaffold_id = liste[0]
            anvio_id = liste[1]
            anvio_updated_id = anvio_id
            refineM_outlier = liste[2]

            if anvio_id == 'Unbinned' :
                refineM_outlier = 'NULL'

            if liste[3] == 'Na' :
                anvio_length = 'NULL'
            else:
                anvio_length = int(liste[3])
            
            if liste[4] == 'Na' :
                anvio_gc = 'NULL'
            else:
                anvio_gc = float(liste[4])
                
            if liste[5] == 'Na' :
                anvio_nb_splits = 'NULL'
            else:
                anvio_nb_splits = int(liste[5])

            if liste[6] == 'Na' :
                anvio_coverage = 'NULL'
            else:
                anvio_coverage = float(liste[6])

            if liste[7] == 'Na':
                anvio_taxonomy = 'unclassified'
            else:
                anvio_taxonomy = liste[7]

            refineM_length = int(liste[8])
            refineM_coverage = float(liste[9])
            refineM_gc = float(liste[10])
            refineM_nb_genes = int(liste[11])
            refineM_coding_length = int(liste[12])
            refineM_domain = liste[13]
            refineM_phylum = liste[18]
            refineM_class = liste[23]
            refineM_order = liste[28]
            refineM_family = liste[33]
            refineM_genus = liste[38]
            refineM_species = liste[43]
            
            if anvio_id == 'Unbinned' and refineM_length < 3000 :
                continue

            row_id = create_scaffold(conn,[ scaffold_id , anvio_id , anvio_updated_id , anvio_gc , anvio_nb_splits , anvio_coverage , anvio_taxonomy , anvio_length , refineM_outlier , refineM_length , refineM_coverage , refineM_gc , refineM_nb_genes , refineM_coding_length , refineM_domain , refineM_phylum , refineM_class , refineM_order , refineM_family , refineM_genus , refineM_species ])
            #print(row_id)

        file.close()




def main():
    sample = 'Esp5_M_AM1'
    project = 'BSI'

    directory = '/env/cns/proj/projet_CSD/scratch/assemblies/'+sample
    database = directory+'/'+'refinedBins'+'/'+'REFINEDBINS.db'
    bins_summary_filename = directory+'/'+'refinedBins'+'/'+'output'+'/'+'Bins_summary.tsv'
    directory = directory+'/'+'refinedBins'+'/'+'output'

    sql_create_scaffolds_table = """ CREATE TABLE IF NOT EXISTS scaffolds (
                                    scaffold_id text PRIMARY KEY,
                                    refineM_outlier text NOT NULL,
                                    refineM_length integer NOT NULL,
                                    refineM_coverage real NOT NULL,
                                    refineM_gc real NOT NULL,
                                    refineM_coding_length integer NOT NULL,
                                    refineM_nb_genes integer NOT NULL,
                                    refineM_domain text,
                                    refineM_phylum text,
                                    refineM_order text,
                                    refineM_class text,
                                    refineM_family text,
                                    refineM_genus text,
                                    refineM_species text,
                                    anvio_updated_id text NOT NULL,
                                    anvio_id text NOT NULL,
                                    anvio_gc real,
                                    anvio_nb_splits integer,
                                    anvio_coverage real,
                                    anvio_taxonomy text,
                                    anvio_length integer,
                                    FOREIGN KEY (scaffold_id) REFERENCES bins (anvio_id)
                                ); """


    sql_create_bins_table = """ CREATE TABLE IF NOT EXISTS bins (
                                        anvio_id text PRIMARY KEY,
                                        name text,
                                        anvio_coverage integer,
                                        anvio_length integer,
                                        anvio_contig_nb integer,
                                        anvio_gc real,
                                        anvio_N50 integer,
                                        anvio_completeness real,
                                        anvio_contamination real,
                                        anvio_taxonomy text,
                                        gtdb_taxonomy text,
                                        gtdb_fastani_reference text,
                                        gtdb_classification_method text,
            	                        checkM_nb_predicted_genes integer,
                                        checkM_translation_table integer,
	                                checkM_coding_density real,
	                                checkM_completeness real,
	                                checkM_contamination real,
                                        anvio_bin integer NOT NULL
                                    ); """


    # create a database connection
    conn = create_connection(database)

    # create tables
    if conn is not None:
        # create bins table
        print('create bins table...')
        create_table(conn, sql_create_bins_table)
        populate_bins_table(conn , bins_summary_filename , sample , project)
        print('done')
        # create scaffolds table
        print('create scaffolds table...')
        create_table(conn, sql_create_scaffolds_table)
        populate_scaffolds_table(conn , directory)
        print('done')
    else:
        print("Error! cannot create the database connection.")


if __name__ == '__main__':
    main()

