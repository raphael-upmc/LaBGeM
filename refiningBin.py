#! /env/cns/proj/agc/scratch/conda/miniconda3/bin/python

import os,sys,re
from collections import defaultdict
from Bio import SeqIO
import argparse
import ast
import json
from xlsxwriter.workbook import Workbook
import csv
import sqlite3
from sqlite3 import Error


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



def runningRefineM(refiningBins_directory,refineM_dir, bin_dir,contig_filename,bam_filename,cpu) :

    print('refineM...')
    os.mkdir(refineM_dir)

    # Removing contaminations based on genomic properties
    genomic_dir = refineM_dir+'/'+'genomicProperties'
    if os.path.exists(genomic_dir) :
        sys.exit(genomic_dir+' already exist, please remove it first')
    os.mkdir(genomic_dir)

    stat_dir = genomic_dir+'/'+'stats'
    if os.path.exists(stat_dir) :
        sys.exit(stat_dir+' already exist, please remove it first')
    os.mkdir(stat_dir)


    cmd = 'source activate refineM-0.1.2 && refinem scaffold_stats -c '+cpu+' '+contig_filename+' '+bin_dir+' '+stat_dir+' '+bam_filename+' >/dev/null 2>&1'
    print(cmd)
    status = os.system(cmd)
    print('status: '+str(status))
    if not os.path.exists(stat_dir+'/scaffold_stats.tsv'):
        sys.exit('something went wrong with refinem scaffold_stats, exit')


    outliers_dir = genomic_dir+'/'+'outliers'
    if os.path.exists(outliers_dir) :
        sys.exit(outliers_dir+' already exist, please remove it first')
    os.mkdir(outliers_dir)

    cmd = 'source activate refineM-0.1.2 && refinem outliers '+stat_dir+'/scaffold_stats.tsv'+' '+outliers_dir+' >/dev/null 2>&1'
    print(cmd)
    status = os.system(cmd)
    print('status: '+str(status))
    if not status == 0:
        sys.exit('something went wrong with refinem scaffold_stats, exit')


    filter_dir = genomic_dir+'/'+'filter'
    if os.path.exists(filter_dir) :
        sys.exit(filter_dir+' already exist, please remove it first')
    os.mkdir(filter_dir)

    genomicOutliers_filename = outliers_dir+'/outliers.tsv'
    cmd = 'source activate refineM-0.1.2 && refinem filter_bins '+bin_dir+' '+genomicOutliers_filename+' '+filter_dir+' >/dev/null 2>&1'
    print(cmd)
    status = os.system(cmd)
    print('status: '+str(status))
    if not status == 0:
        sys.exit('something went wrong with refinem filter_bins, exit')
    print('\n\n')


    # Removing contaminations based on the contig taxonomy

    gene_dir = anvio_directory+'/'+'genes'
    if os.path.exists(gene_dir) :
        sys.exit(gene_dir_dir+' already exist, please remove it first')
    os.mkdir(gene_dir)

    # replace the call cmd line by an housemade function to get the cds,gff,protein files
    cmd = 'source activate refineM-0.1.2 && refinem call_genes -c '+cpu+' '+bin_dir+' '+gene_dir+' >/dev/null 2>&1'
    print(cmd)
    status = os.system(cmd)
    print('status: '+str(status))
    if not status == 0:
        sys.exit('something went wrong with refinem call_genes, exit')

    reference_taxonomy_filename = '/env/cns/proj/agc/scratch/proj/GTDB/gtdb_r95_taxonomy.2020-07-30.tsv'
    reference_db_filename = '/env/cns/proj/agc/scratch/proj/GTDB/gtdb_r95_protein_db.2020-07-30.faa.dmnd'

    taxo_dir = refineM_dir+'/'+'taxonomy'
    if os.path.exists(taxo_dir) :
        sys.exit(taxo_dir+' already exist, please remove it first')
    os.mkdir(taxo_dir)

    taxoProfile_dir = taxo_dir+'/'+'profiles'
    if os.path.exists(taxoProfile_dir) :
        sys.exit(taxoProfile_dir+' already exist, please remove it first')
    os.mkdir(taxoProfile_dir)

    log_filename = refineM_dir+'/'+'taxonomy'+'/'+'taxon_profile.log'
    cmd = 'source activate refineM-0.1.2 && refinem taxon_profile -c '+cpu+' '+gene_dir+' '+stat_dir+'/scaffold_stats.tsv'+' '+reference_db_filename+' '+reference_taxonomy_filename+' '+taxoProfile_dir+' >'+log_filename+' 2>&1'
    print(cmd)
    status = os.system(cmd)
    print('status: '+str(status))
    if not status == 0:
        sys.exit('something went wrong with refinem taxon_profile, exit')

    taxoOutlier_dir = taxo_dir+'/'+'outliers'
    if os.path.exists(taxoOutlier_dir) :
        sys.exit(taxoOutlier_dir+' already exist, please remove it first')
    os.mkdir(taxoOutlier_dir)


    taxoOutliers_filename = taxoOutlier_dir+'/'+'taxon_filter.tsv'
    cmd = 'source activate refineM-0.1.2 && refinem taxon_filter -c '+cpu+' '+taxoProfile_dir+' '+taxoOutliers_filename+' >/dev/null 2>&1'
    print(cmd)
    status = os.system(cmd)
    print('status: '+str(status))
    if not status == 0:
        sys.exit('something went wrong with refinem taxon_filter, exit')

    taxoFilter_dir = taxo_dir+'/'+'filters'
    if os.path.exists(taxoFilter_dir) :
        sys.exit(taxoFilter_dir+' already exist, please remove it first')
    os.mkdir(taxoFilter_dir)

    cmd = 'source activate refineM-0.1.2 && refinem filter_bins '+bin_dir+' '+taxoOutliers_filename+' '+taxoFilter_dir+' >/dev/null 2>&1'
    print(cmd)
    status = os.system(cmd)
    print('status: '+str(status))
    if not status == 0:
        sys.exit('something went wrong with refinem filter_bins, exit')



def sortingBAM(bam_filename,bai_filename) :
    tmp_bam_filename = bam_filename+'.unsorted'
    os.rename(bam_filename,tmp_bam_filename)
    cmd = 'module load samtools/1.10.2 && samtools sort -@ '+str(cpu)+' -o '+bam_filename+' -O BAM '+tmp_bam_filename+' >/dev/null 2>&1'
    print(cmd)
    status = os.system(cmd)
    print('status: '+str(status))

    # creating the index file
    cmd = 'module load samtools/1.10.2 && samtools index '+bam_filename+' >/dev/null 2>&1'
    print(cmd)
    status = os.system(cmd)
    print('status: '+str(status)+'\n')
    if not status == 0 :
        sys.exit('something went wrong with samtools index, exit.')

    os.remove(tmp_bam_filename)
    if os.path.exists(tmp_bam_filename) :
        sys.exit('something went wrong with os.remove, exit')

    if not os.path.exists(bai_filename):
        sys.exit('something went wrong with samtools sort, exit')


def getProjectSampleNames(scaffold_filename) :
    file = open(scaffold_filename,'r')
    for line in file :
        if re.match('>',line) :
            defline = line.rstrip().replace('>','')
            project,sample,contig = defline.split('__')
            break
        else:
            continue
    file.close()
    return project,sample


def update(directory) :

    scaffold_filename =  directory+'/'+'assembly'+'/'+'megahit.contigs.renamed.fa'
    protein_filename =   directory+'/'+'assembly'+'/'+'proteins.anvio.tab'
    contigDb_filename =   directory+'/'+'assembly'+'/'+'contigs.db'
    bam_filename = directory+'/'+'assembly'+'/'+'bt2'+'/'+'megahit.contigs.renamed.fa.bam'
    bai_filename = directory+'/'+'assembly'+'/'+'bt2'+'/'+'megahit.contigs.renamed.fa.bam.bai'
    profileDb_filename = directory+'/'+'assembly'+'/'+'anvio'+'/'+'PROFILE.db'

    project,sample = getProjectSampleNames(scaffold_filename)


    for root, dirs, files in os.walk(directory+'/'+'assembly'+'/'+'bt2', topdown = False):
        for filename in files :
            if re.search(r'sorted',filename) :
                if re.search(r'.sorted.fake.bam$',filename) :
                    continue
                if re.search(r'.sorted.fake.bam.bai$',filename) :
                    continue

                if re.search(r'.sorted.bam$',filename) :
                    anvio_bam_filename = root+'/'+filename
#                    print(anvio_bam_filename)

                if re.search(r'.sorted.bam.bai$',filename) :
                    anvio_bai_filename = root+'/'+filename
#                    print(anvio_bai_filename)

#    print('\n\n')
#    print('looking for anvio contig filename...')
    for root, dirs, files in os.walk(directory+'/'+'assembly', topdown = True):
        if os.path.abspath(root) != os.path.abspath(directory+'/'+'assembly') :
            continue
        for filename in files :
            if re.search(r'^megahit.contigs.renamed',filename) and filename != 'megahit.contigs.renamed.fa' :
                anvio_contig_filename = root+'/'+filename
#    print(anvio_contig_filename)
#    sys.exit()

    datatable_dir = directory+'/'+'assembly'+'/'+'datatables'
    # if datatables isn't prensent, create it #
    if not os.path.exists(datatable_dir) :
        os.mkdir(datatable_dir)


    taxo_anvio_filename = datatable_dir+'/'+'taxon_names.txt'
    if not os.path.exists(taxo_anvio_filename) :
        cmd = 'source activate anvio-6.2 && anvi-export-table '+contigDb_filename+' --table taxon_names -o '+taxo_anvio_filename+' >/dev/null 2>&1'
        print(cmd)
        status = os.system(cmd)
        print('status: '+str(status)+'\n')
        if not status == 0 :
            sys.exit('something went wrong with anvi-export-table, exit.')


    gene_taxo_anvio_filename = datatable_dir+'/'+'genes_taxonomy.txt'
    if not os.path.exists(gene_taxo_anvio_filename) :
        cmd = 'source activate anvio-6.2 && anvi-export-table '+contigDb_filename+' --table genes_taxonomy -o '+gene_taxo_anvio_filename+' >/dev/null 2>&1'
        print(cmd)
        status = os.system(cmd)
        print('status: '+str(status)+'\n')
        if not status == 0 :
            sys.exit('something went wrong with anvi-export-table, exit.')


    basic_info_contigs_filename = datatable_dir+'/'+'contigs_basic_info.txt'
    if not os.path.exists(basic_info_contigs_filename) :
        cmd = 'source activate anvio-6.2 && anvi-export-table '+contigDb_filename+' --table contigs_basic_info -o '+basic_info_contigs_filename+' >/dev/null 2>&1'
        print(cmd)
        status = os.system(cmd)
        print('status: '+str(status)+'\n')
        if not status == 0 :
            sys.exit('something went wrong with anvi-export-table, exit.')


    coverage_contigs_filename = datatable_dir+'/'+'contigs_coverage_info.txt'
    if not os.path.exists(coverage_contigs_filename) :
        cmd = 'source activate anvio-6.2 && anvi-export-splits-and-coverages -p '+profileDb_filename+' -c '+contigDb_filename+' -o '+datatable_dir+' -O '+'tmp'+' --report-contigs'+' >/dev/null 2>&1'
        print(cmd)
        status = os.system(cmd)
        print('status: '+str(status)+'\n')
        if not status == 0 :
            sys.exit('something went wrong with anvi-export-splits-and-coverages, exit.')
        os.rename(datatable_dir+'/'+'tmp-COVs.txt',coverage_contigs_filename)
        os.remove(datatable_dir+'/'+'tmp-CONTIGS.fa' )

    # if the bam_filename isn't sorted, do #
    if not os.path.exists(bai_filename) :
        sortingBAM(bam_filename,bai_filename)

    # creating a json file #
    config_filename = directory+'/'+'assembly'+'/'+'info.json'
    output = open(config_filename,'w')
    output.write('{'+'\n')
    output.write('\t\"project\":\"'+project+'\",\n')
    output.write('\t\"sample\":\"'+sample+'\",\n')
    output.write('\t\"directory\":\"'+directory+'\",\n')
    output.write('\t\"assembly_cmd_line\":\"'+'Na'+'\",\n') # ??
    output.write('\t\"assembly_directory\":\"'+directory+'/'+'assembly'+'\",\n')
    output.write('\t\"assembly_contig_filename\":\"'+directory+'/'+'assembly'+'/'+'megahit.contigs.renamed.fa'+'\",\n')
    output.write('\t\"assembly_bam_filename\":\"'+bam_filename+'\",\n')
    output.write('\t\"assembly_bai_filename\":\"'+bai_filename+'\",\n')
    output.write('\t\"assembly_protein_filename\":\"'+directory+'/'+'assembly'+'/'+'proteins.faa'+'\",\n')
    output.write('\t\"assembly_gene_filename\":\"'+directory+'/'+'assembly'+'/'+'genes.fna'+'\",\n')
    output.write('\t\"anvio_bai_filename\":\"'+anvio_bai_filename+'\",\n')
    output.write('\t\"anvio_bam_filename\":\"'+anvio_bam_filename+'\",\n')
    output.write('\t\"anvio_contig_filename\":\"'+directory+'/'+'assembly'+'/'+'megahit.contigs.renamed.fa'+'\",\n') # ??
    output.write('\t\"anvio_protein_filename\":\"'+directory+'/'+'assembly'+'/'+'proteins.anvio.tab'+'\",\n')
    output.write('\t\"anvio_contigDb_filename\":\"'+directory+'/'+'assembly'+'/'+'contigs.db'+'\",\n')
    output.write('\t\"anvio_profileDb_filename\":\"'+directory+'/'+'assembly'+'/'+'anvio'+'/'+'PROFILE.db'+'\",\n')
    output.write('\t\"anvio_coverage_contigs_filename\":\"'+coverage_contigs_filename+'\",\n')
    output.write('\t\"anvio_basic_info_contigs_filename\":\"'+basic_info_contigs_filename+'\",\n')
    output.write('\t\"anvio_gene_taxo_anvio_filename\":\"'+gene_taxo_anvio_filename+'\",\n')
    output.write('\t\"anvio_taxo_filename\":\"'+taxo_anvio_filename+'\"\n')
    output.write('}'+'\n')
    output.close()



def runningGTDBtk(gtdbtk_dir,bin_dir,cpu) :
    print('GTDB-tk...')
    os.mkdir(gtdbtk_dir)
    os.mkdir(gtdbtk_dir+'/'+'bins')
    os.mkdir(gtdbtk_dir+'/'+'output')
    
    for root, dirs, files in os.walk(bin_dir, topdown = False):
        for filename in files :
            if re.search(r'.fna',filename) :
                print(filename)
                binName = filename.replace('.fna','')
                if binName == 'Unbinned' :
                    continue
                if re.match(r'Euk',binName) :
                    continue
                print(root+'/'+filename)
                os.symlink(root+'/'+filename,gtdbtk_dir+'/'+'bins'+'/'+binName+'.fna')
    cmd = 'source activate gtdbtk-1.4.0 && gtdbtk classify_wf  --cpus '+cpu+' --genome_dir '+gtdbtk_dir+'/'+'bins'+' --out_dir '+gtdbtk_dir+'/'+'output > '+gtdbtk_dir+'/'+'gtdbtk.log 2>&1'
    print(cmd)
    status = os.system(cmd)
    print('status: '+str(status)+'\n')
    if not status == 0 :
        sys.exit('something went wrong with gtdbtk classify_wf, exit.')



def runningCheckM(checkm_dir,bin_dir,cpu) :
    print('checkM...')
    os.mkdir(checkm_dir)
    os.mkdir(checkm_dir+'/'+'bins')
    os.mkdir(checkm_dir+'/'+'output')
    
    for root, dirs, files in os.walk(bin_dir, topdown = False):
        for filename in files :
            if re.search(r'.fna',filename) :
                print(filename)
                binName = filename.replace('.fna','')
                if binName == 'Unbinned' :
                    continue
                if re.match(r'Euk',binName) :
                    continue
                print(root+'/'+filename)
                os.symlink(root+'/'+filename,checkm_dir+'/'+'bins'+'/'+binName+'.fna')
    cmd = 'source activate checkM-1.1.3 && checkm lineage_wf -t '+cpu+' '+checkm_dir+'/'+'bins'+' '+checkm_dir+'/'+'output > '+checkm_dir+'/'+'checkm.log 2>&1'
    print(cmd)
    status = os.system(cmd)
    print('status: '+str(status)+'\n')
    if not status == 0 :
        sys.exit('something went wrong with checkm lineage_wf, exit.')


def writingOutput(json_data, refiningBins_directory , anvio_scaffold2info, anvio_scaffold2taxonomy , partialScaffold2bin) :

    #################################
    # assembly and collection info #
    #################################

    #######################
    # anvio sample summary
    collection_filename = refiningBins_directory+'/ANVIO/collection.txt'
    file = open(collection_filename,'r')
    for line in file :
        collection = line.rstrip()
    file.close()

    bin2anvio_summary = dict()
    anvio_summary_filename = refiningBins_directory+'/ANVIO/SAMPLES-SUMMARY/bins_summary.txt'
    file = open(anvio_summary_filename,'r')
    headerList = next(file).rstrip().split('\t')
    print(headerList)
    for line in file :
        line = line.rstrip()
        liste = line.split('\t')
        print(liste)
        binName = liste[0]
        if binName not in bin2anvio_summary :
            bin2anvio_summary[ binName ] = dict()

        for i in range(1,len(headerList)) :
            header = 'Anvio_'+headerList[i].replace(' ','_')
            try :
                feature = liste[i]
            except:
                feature = 'Na'
            bin2anvio_summary[binName][header] = feature
        print(bin2anvio_summary[binName])
    file.close()

    anvio_genome_coverage_filename = refiningBins_directory+'/ANVIO/SAMPLES-SUMMARY/bins_across_samples/mean_coverage.txt'
    file = open(anvio_genome_coverage_filename,'r')
    header = next(file).rstrip().split('\t')
    for line in file :
        line = line.rstrip()
        liste = line.split('\t')
        binName = liste[0]
        coverage = liste[1]
        bin2anvio_summary[binName]['Anvio_mean_coverage'] = coverage
    file.close()


    #########
    # checkM

    checkmHeaderSet = set()
    bin2checkm = dict()
    checkm_filename = refiningBins_directory+'/CheckM/output/storage/bin_stats_ext.tsv'
    file = open(checkm_filename,'r')
    for line in file :
        line = line.rstrip()
        liste = line.split('\t')
        bin2checkm[ liste[0] ] = ast.literal_eval(liste[1])
        for key,value in bin2checkm[ liste[0] ].items() :
            checkmHeaderSet.add(key)
    file.close()


    #######
    # gtdb

    bin2gtdb = dict()
    gtdb_arc_filename = refiningBins_directory+'/GTDB-tk/output/gtdbtk.ar122.summary.tsv'
    gtdb_bac_filename = refiningBins_directory+'/GTDB-tk/output/gtdbtk.bac120.summary.tsv'

    for filename in [gtdb_bac_filename,gtdb_arc_filename] :
#        print(filename)
        if not os.path.exists(filename) :
            continue
        else:
            file = open(filename,'r')
            gtdb_headerList = next(file).rstrip().split('\t')
            for line in file :
                line = line.rstrip()
                liste = line.split('\t')
                binName = liste[0]
                if binName not in bin2gtdb :
                    bin2gtdb[ binName ] = dict()

                for i in range(0,len(gtdb_headerList)) :
                    header = gtdb_headerList[i]
                    feature = liste[i]
                    #print(str(header)+'\t'+str(feature))
                    bin2gtdb[binName][header] = feature
            file.close()
#    print(bin2gtdb)



    ###############
    # collections #
    ###############

    filenameList = list()
    output_dir = refiningBins_directory+'/'+'output' 
    os.mkdir(output_dir)

    print('writting outputs in'+output_dir+'...')
    filenameList.append(output_dir+'/'+'Sample_summary.tsv')
    filenameList.append(output_dir+'/'+'Bins_summary.tsv')
    filenameList.append(output_dir+'/'+'Anvio_summary.tsv')
    filenameList.append(output_dir+'/'+'CheckM.tsv')
    filenameList.append(output_dir+'/'+'GTDBtk.tsv')


    # Anvio-summary
    output = open(output_dir+'/'+'Anvio_summary.tsv','w')
    headerAnvioList = ['Anvio_taxon','Anvio_mean_coverage','Anvio_total_length','Anvio_num_contigs','Anvio_N50','Anvio_GC_content','Anvio_percent_completion','Anvio_percent_redundancy']
    output.write('Bin'+'\t'+'\t'.join(headerAnvioList)+'\n')
    for binName,key2feature in bin2anvio_summary.items() :
        output.write(binName)
        print(binName)
        for header in headerAnvioList :
            print(header+'\t'+key2feature[header])
            output.write('\t'+key2feature[header])
        output.write('\n')
    output.close()
    
    # CheckM
    output = open(output_dir+'/'+'CheckM.tsv','w')
    output.write( 'Bin'+'\t'+'\t'.join(sorted(list(checkmHeaderSet)))+'\n')
    for binName,key2feature in bin2checkm.items() :
        output.write(binName)
        for header in sorted(list(checkmHeaderSet)) :
            if header in bin2checkm[binName] :
                output.write('\t'+str(bin2checkm[binName][header]))
            else:
                output.write('\t'+'Na')
        output.write('\n')
    output.close()

    # GTDB
    print(gtdb_headerList)
    output = open(output_dir+'/'+'GTDBtk.tsv','w')
    output.write( 'Bin'+'\t'+'\t'.join(gtdb_headerList)+'\n' )
    for binName,key2feature in bin2gtdb.items() :
        output.write(binName)
        for header in gtdb_headerList :
            if header in bin2gtdb[binName] :
                output.write('\t'+bin2gtdb[binName][header])
            else:
                output.write('\t'+'Na')                                
        output.write('\n')
    output.close()

    # Collection
    output = open(output_dir+'/'+'Sample_summary.tsv','w')
    output.write('Project: '+'\t'+json_data['project']+'\n')
    output.write('Sample: '+'\t'+json_data['sample']+'\n')
    output.write('Collection: '+'\t'+collection+'\n')
    output.write('Assembly directory: '+'\t'+json_data['assembly_directory']+'\n')
    output.write('Working directory: '+'\t'+refiningBins_directory+'\n')
    output.close()

    output = open(output_dir+'/'+'Bins_summary.tsv','w')
    output.write('Bin'+'\t'+'\t'.join(headerAnvioList)+'\t'+'\t'.join(['CheckM_#_predicted_genes','CheckM_Translation_table','CheckM_Coding_density','CheckM_Completeness','CheckM_Contamination'])+'\t'+'\t'.join( ['Gtdb_classification','Gtdb_fastani_reference','Gtdb_classification_method'] )+'\n')
    for binName,key2feature in bin2anvio_summary.items() :
        output.write(binName)
        for header in headerAnvioList :
            output.write('\t'+key2feature[header])

        for header in ['# predicted genes','Translation table','Coding density','Completeness','Contamination'] :
            if binName in bin2checkm and header in bin2checkm[binName] :
                output.write('\t'+str(bin2checkm[binName][header]))
            else :
                output.write('\t'+'Na')

        for header in ['classification','fastani_reference','classification_method'] :
            if binName in bin2gtdb and header in bin2gtdb[binName] :
                output.write('\t'+str(bin2gtdb[binName][header]))
            else :
                output.write('\t'+'Na')
        output.write('\n')
    output.close()

    scaffold_filename = output_dir+'/'+'partial_scaffolds.tsv'
    filenameList.append(scaffold_filename)
    output = open(scaffold_filename,'w')
    output.write('scaffold'+'\t'+'associated_bins'+'\t'+'final_bin'+'\n')
    for scaffold,liste in partialScaffold2bin.items() :
        if len(liste) == 1 :
            output.write(scaffold+'\t'+','.join(list(liste))+'\t'+list(liste)[0]+'\n')
        else:
            output.write(scaffold+'\t'+','.join(list(liste))+'\t'+'Unbinned'+'\n')
    output.close()


    ##########################
    # creating the bin files #
    ##########################

    ###########
    # refineM
    genomicOutliers_filename = refiningBins_directory+'/refineM/genomicProperties/outliers/outliers.tsv'
    taxoOutliers_filename = refiningBins_directory+'/refineM/taxonomy/outliers/taxon_filter.tsv'
    taxoProfile_dir = refiningBins_directory+'/refineM/taxonomy/profiles/bin_reports'
    genomicStat_filename =  refiningBins_directory+'/refineM'+'/'+'genomicProperties/stats/scaffold_stats.tsv'

    print('refineM...')
    outliersSet = set()
    scaffold2outlier = dict()
    file = open(genomicOutliers_filename,'r')
    header = next(file)
    for line in file :
        line = line.rstrip()
        liste = line.split()
        scaffold = liste[0]
        outliersSet.add(scaffold)
        scaffold2outlier[scaffold] = liste[3]
    file.close()


    scaffold2genomicProperties = dict()
    file = open(genomicStat_filename,'r')
    header = next(file)
    for line in file :
        line = line.rstrip()
        liste = line.split('\t')
        scaffold = liste[0]
        gc = str( float(liste[2]) / 100.0 )
        length = liste[3]
        coverage = liste[4]
        scaffold2genomicProperties[scaffold] = [length,coverage,gc]
    file.close()


    outliersTaxoSet = set()
    file = open(taxoOutliers_filename,'r')
    for line in file :
        line = line.rstrip()
        if re.match(r'#',line) :
            continue
        liste = line.split()
        scaffold = liste[0]
        outliersTaxoSet.add(scaffold)
        if scaffold in scaffold2outlier :
            scaffold2outlier[scaffold] += ',TAXO'
        else:
            scaffold2outlier[scaffold] = 'TAXO'
    file.close()


    for scaffold,liste in partialScaffold2bin.items() :
        if scaffold in scaffold2outlier :
            scaffold2outlier[scaffold] += ',PARTIAL ('+','.join(list(liste))+')'
        else:
            scaffold2outlier[scaffold] = 'PARTIAL ('+','.join(list(liste))+')'

    refineM_scaffold2info = dict()
    for root, dirs, files in os.walk(taxoProfile_dir, topdown = False):
        for filename in files :
            if re.search(r'.scaffolds.tsv',filename) :
                binName = filename.replace('.scaffolds.tsv','')
                print(binName+'\t'+root+'/'+filename)
                file = open(root+'/'+filename,'r')
                header = next(file).rstrip()
                for line in file :
                    line = line.rstrip()
                    liste = line.split('\t')
                    scaffold = liste[0]
                    refineM_scaffold2info[ scaffold ] = '\t'.join( line.split('\t')[5:] )
                file.close()


    for scaffold,liste in partialScaffold2bin.items() :
        scaffold2outlier[scaffold] += ',PARTIAL ('+','.join(list(liste))+')'



    ############
    # bin files

    scaffold2length = dict()
    bin2scaffold = defaultdict(set)
    bin_directory =  directory+'/'+'refinedBins'+'/'+'ANVIO'+'/'+'bins'
    for root, dirs, files in os.walk(bin_directory, topdown=False):
        for filename in files :
            binName = filename.replace('.fna','')
            for record in SeqIO.parse(root+'/'+filename,'fasta') :
                scaffold2length[record.id] = len(record)
                bin2scaffold[binName].add(record.id)



    headerRefineM_taxo = '\t'.join( header.split('\t')[5:] )

    for binName,scaffoldSet in bin2scaffold.items() :
        output_filename = output_dir+'/'+binName+'.tsv'
        filenameList.append(output_filename)

        print('\n'+output_filename+'\t'+str(len(scaffoldSet)))

        output = open(output_filename,'w')
        output.write('scaffold'+'\t'+'bin'+'\t'+'refineM_outlier'+'\t'+'anvio_length'+'\t'+'anvio_gc'+'\t'+'anvio_nb_splits'+'\t'+'anvio_coverage'+'\t'+'anvio_taxonomy'+'\t'+'refineM_length'+'\t'+'refineM_coverage'+'\t'+'refineM_gc'+'\t'+headerRefineM_taxo+'\n')
        for scaffold,length in sorted(scaffold2length.items() , key = lambda x:x[1] , reverse = True ) :

            if scaffold not in scaffoldSet :
                continue

            if scaffold in anvio_scaffold2taxonomy :
                taxonomy_anvio = anvio_scaffold2taxonomy[scaffold]
            else:
                taxonomy_anvio = 'Na'
                
            if scaffold in anvio_scaffold2info :
                info_anvio = '\t'.join(anvio_scaffold2info[ scaffold ])
            else:
                info_anvio = 'Na\tNa\tNa\tNa'

            if scaffold in scaffold2genomicProperties :
                info_refineM = '\t'.join(scaffold2genomicProperties[ scaffold ])
            else:
                info_refineM = 'Na\tNa\tNa'

            if scaffold in refineM_scaffold2info :
                taxo_refineM = refineM_scaffold2info[scaffold]
            else:
                taxo_refineM = '0'+'\t'+'0'+'\t'+'unclassified'+'\t'+'na'+'\t'+'na'+'\t'+'na'+'\t'+'na'+'\t'+'unclassified'+'\t'+'na'+'\t'+'na'+'\t'+'na'+'\t'+'na'+'\t'+'unclassified'+'\t'+'na'+'\t'+'na'+'\t'+'na'+'\t'+'na'+'\t'+'unclassified'+'\t'+'na'+'\t'+'na'+'\t'+'na'+'\t'+'na'+'\t'+'unclassified'+'\t'+'na'+'\t'+'na'+'\t'+'na'+'\t'+'na'+'\t'+'unclassified'+'\t'+'na'+'\t'+'na'+'\t'+'na'+'\t'+'na'+'\t'+'unclassified'+'\t'+'na'+'\t'+'na'+'\t'+'na'+'\t'+'na'

            if scaffold in scaffold2outlier :
                outlier = scaffold2outlier[scaffold]
            else :
                outlier = '-'

            output.write(scaffold+'\t'+binName+'\t'+outlier+'\t'+info_anvio+'\t'+taxonomy_anvio+'\t'+info_refineM+'\t'+taxo_refineM+'\n')
        output.close()



    ###########################
    # building the excel file #
    ###########################

    # https://towardsdatascience.com/writing-to-excel-with-python-micropython-42cf9541c101
    xlsx_filename = output_dir+'/'+'Collection.xlsx'

    # Create an XlsxWriter workbook object and add a worksheet.
    book = Workbook(xlsx_filename)

    #print(filenameList)
    for filename in filenameList :
        basename = os.path.basename(filename).replace('.tsv','')
        sheet = book.add_worksheet(basename)
        # print(filename)

        cpt = 0
        file = open(filename, 'r')
        for line in file: # Read the row data from the TSV file and write it to the XLSX file.                
            line = line.strip()
            liste = line.split('\t')
            sheet.write_row(cpt, 0, liste)
            cpt += 1
        file.close()

    # Close the XLSX file.
    book.close()



def writingOutputProkaryote(json_data, refiningBins_directory , anvio_scaffold2info, anvio_scaffold2taxonomy , partialScaffold2bin) : # if no prok in the metagenome

    #################################
    # assembly and collection info #
    #################################

    #######################
    # anvio sample summary
    collection_filename = refiningBins_directory+'/ANVIO/collection.txt'
    file = open(collection_filename,'r')
    for line in file :
        collection = line.rstrip()
    file.close()

    bin2anvio_summary = dict()
    anvio_summary_filename = refiningBins_directory+'/ANVIO/SAMPLES-SUMMARY/bins_summary.txt'
    file = open(anvio_summary_filename,'r')
    headerList = next(file).rstrip().split('\t')
    for line in file :
        line = line.rstrip()
        liste = line.split('\t')
        binName = liste[0]
        if binName not in bin2anvio_summary :
            bin2anvio_summary[ binName ] = dict()

        for i in range(1,len(headerList)) :
            header = 'Anvio_'+headerList[i].replace(' ','_')
            feature = liste[i]
            bin2anvio_summary[binName][header] = feature
    file.close()

    anvio_genome_coverage_filename = refiningBins_directory+'/ANVIO/SAMPLES-SUMMARY/bins_across_samples/mean_coverage.txt'
    file = open(anvio_genome_coverage_filename,'r')
    header = next(file).rstrip().split('\t')
    for line in file :
        line = line.rstrip()
        liste = line.split('\t')
        binName = liste[0]
        coverage = liste[1]
        bin2anvio_summary[binName]['Anvio_mean_coverage'] = coverage
    file.close()


    ###############
    # collections #
    ###############

    filenameList = list()
    output_dir = refiningBins_directory+'/'+'output' 
    os.mkdir(output_dir)

    print('writting outputs in'+output_dir+'...')
    filenameList.append(output_dir+'/'+'Sample_summary.tsv')
    filenameList.append(output_dir+'/'+'Bins_summary.tsv')
    filenameList.append(output_dir+'/'+'Anvio_summary.tsv')

    # Anvio-summary
    output = open(output_dir+'/'+'Anvio_summary.tsv','w')
    headerAnvioList = ['Anvio_taxon','Anvio_mean_coverage','Anvio_total_length','Anvio_num_contigs','Anvio_N50','Anvio_GC_content','Anvio_percent_completion','Anvio_percent_redundancy']
    output.write('Bin'+'\t'+'\t'.join(headerAnvioList)+'\n')
    for binName,key2feature in bin2anvio_summary.items() :
        output.write(binName)
        for header in headerAnvioList :
            output.write('\t'+key2feature[header])
        output.write('\n')
    output.close()
    
    # Collection
    output = open(output_dir+'/'+'Sample_summary.tsv','w')
    output.write('Project: '+'\t'+json_data['project']+'\n')
    output.write('Sample: '+'\t'+json_data['sample']+'\n')
    output.write('Collection: '+'\t'+collection+'\n')
    output.write('Assembly directory: '+'\t'+json_data['assembly_directory']+'\n')
    output.write('Working directory: '+'\t'+refiningBins_directory+'\n')
    output.close()

    output = open(output_dir+'/'+'Bins_summary.tsv','w')
    output.write('Bin'+'\t'+'\t'.join(headerAnvioList)+'\t'+'\t'.join(['CheckM_#_predicted_genes','CheckM_Translation_table','CheckM_Coding_density','CheckM_Completeness','CheckM_Contamination'])+'\t'+'\t'.join( ['Gtdb_classification','Gtdb_fastani_reference','Gtdb_classification_method'] )+'\n')
    for binName,key2feature in bin2anvio_summary.items() :
        output.write(binName)
        for header in headerAnvioList :
            output.write('\t'+key2feature[header])

        for header in ['# predicted genes','Translation table','Coding density','Completeness','Contamination'] :
            output.write('\t'+'Na')

        for header in ['classification','fastani_reference','classification_method'] :
            output.write('\t'+'Na')
        output.write('\n')
    output.close()

    scaffold_filename = output_dir+'/'+'partial_scaffolds.tsv'
    filenameList.append(scaffold_filename)
    output = open(scaffold_filename,'w')
    output.write('scaffold'+'\t'+'associated_bins'+'\t'+'final_bin'+'\n')
    for scaffold,liste in partialScaffold2bin.items() :
        if len(liste) == 1 :
            output.write(scaffold+'\t'+','.join(list(liste))+'\t'+list(liste)[0]+'\n')
        else:
            output.write(scaffold+'\t'+','.join(list(liste))+'\t'+'Unbinned'+'\n')
    output.close()


    ##########################
    # creating the bin files #
    ##########################

    scaffold2outlier = dict()
    for scaffold,liste in partialScaffold2bin.items() :
        if scaffold in scaffold2outlier :
            scaffold2outlier[scaffold] += ',PARTIAL ('+','.join(list(liste))+')'
        else:
            scaffold2outlier[scaffold] = 'PARTIAL ('+','.join(list(liste))+')'


    ############
    # bin files

    bin_directory =  directory+'/'+'refinedBins'+'/'+'ANVIO'+'/'+'bins'
    for root, dirs, files in os.walk(bin_directory, topdown=False):
        for filename in files :
            binName = filename.replace('.fna','')
            scaffold2length = dict()
            for record in SeqIO.parse(root+'/'+filename,'fasta') :
                scaffold2length[record.id] = len(record)
        
            output_filename = output_dir+'/'+binName+'.tsv'
            filenameList.append(output_filename)
            output = open(output_filename,'w')
            output.write('scaffold'+'\t'+'bin'+'\t'+'refineM_outlier'+'\t'+'anvio_length'+'\t'+'anvio_gc'+'\t'+'anvio_nb_splits'+'\t'+'anvio_coverage'+'\t'+'anvio_taxonomy'+'\n')
            for scaffold,length in sorted(scaffold2length.items() , key = lambda x:x[1] , reverse = True ) :
                if scaffold in anvio_scaffold2taxonomy :
                    taxonomy = anvio_scaffold2taxonomy[scaffold]
                else:
                    taxonomy = 'Na'

                if scaffold in anvio_scaffold2info :
                    info = '\t'.join(anvio_scaffold2info[ scaffold ])
                else:
                    info = 'Na\tNa\tNa\tNa'

                if scaffold in scaffold2outlier :
                    outlier = scaffold2outlier[scaffold]
                else :
                    outlier = '-'
                output.write(scaffold+'\t'+binName+'\t'+outlier+'\t'+info+'\t'+taxonomy+'\n')
            output.close()



    ###########################
    # building the excel file #
    ###########################

    # https://towardsdatascience.com/writing-to-excel-with-python-micropython-42cf9541c101
    xlsx_filename = output_dir+'/'+'Collection.xlsx'

    # Create an XlsxWriter workbook object and add a worksheet.
    book = Workbook(xlsx_filename)

    print(filenameList)
    for filename in filenameList :
        basename = os.path.basename(filename).replace('.tsv','')
        sheet = book.add_worksheet(basename)
        # print(filename)

        cpt = 0
        file = open(filename, 'r')
        for line in file: # Read the row data from the TSV file and write it to the XLSX file.                
            line = line.strip()
            liste = line.split('\t')
            sheet.write_row(cpt, 0, liste)
            cpt += 1
        file.close()

    # Close the XLSX file.
    book.close()




def gettingContigInfo(basic_info_contigs_filename, coverage_contigs_filename ) : # minimum percentage of genes to assign a scaffold to a taxonomic group (default: 20.0)
    scaffold2info  = defaultdict()
    file = open( basic_info_contigs_filename, 'r' )
    header = next(file)
    for line in file :
        line = line.rstrip()
        scaffold,length,gc,nb_splits = line.split('\t')
        scaffold2info[ scaffold ] = [length,gc,nb_splits]
    file.close()

    file = open( coverage_contigs_filename, 'r' )
    header = next(file)
    for line in file :
        line = line.rstrip()
        liste = line.split('\t')
        scaffold = liste[0]
        coverage = liste[1]
        scaffold2info[ scaffold ].append(coverage)
    file.close()
    
    return scaffold2info


def detectingContigTaxonomy(gene_taxo_anvio_filename , taxo_anvio_filename , protein_filename ) : # minimum percentage of genes to assign a scaffold to a taxonomic group (default: 20.0)

    scaffold2genes = defaultdict(set)
    gene2scaffold = dict()
    file = open( protein_filename, 'r' )
    header = next(file)
    for line in file :
        line = line.rstrip()
        liste = line.split('\t')

        geneId = liste[0]
        scaffold = liste[1]
        gene2scaffold[ geneId ] = scaffold
        scaffold2genes[ scaffold ].add(geneId)
    file.close()

    taxoId2taxon = dict()
    file = open( taxo_anvio_filename, 'r' )
    header = next(file)
    for line in file :
        line = line.rstrip()
        liste = line.split('\t')
        taxoId = liste[0]
        del liste[0]
        taxoId2taxon[taxoId] = liste
    file.close()

    geneId2taxoId = dict()
    file = open( gene_taxo_anvio_filename, 'r' )
    header = next(file)
    for line in file :
        line = line.rstrip()
        geneId,taxoId = line.split('\t')
        geneId2taxoId[ geneId ] = taxoId
    file.close()

    scaffold2taxonomy = dict()
    for scaffold,geneSet in scaffold2genes.items() :
        #print()
        #print(scaffold)
        taxo2pct = defaultdict(float)
        for geneId in sorted(geneSet) :
            if geneId in geneId2taxoId :
                taxoId = geneId2taxoId[ geneId ]
                taxon = taxoId2taxon[taxoId][-2]
            else:
                taxon = 'Unknown'

            #print('\t'+geneId+'\t'+str(taxon) )
            taxo2pct[ taxon ] += 1 / float(len(geneSet))
        #print('\t'+str(taxo2pct))
        #print()

        TAXON = ''
        for taxon,pct in sorted(taxo2pct.items(),key=lambda x:x[1], reverse=True) :
            #print('\t\t'+taxon+'\t'+str(pct))
            if pct > 0.2 and taxon != 'Unknown' :
                if TAXON == '' :
                    TAXON = taxon

        if TAXON != '' :
            #print('\t\ttaxon: '+TAXON)
            scaffold2taxonomy[ scaffold ] = TAXON

    return scaffold2taxonomy


# def writingDatabase() :
#     # toto


if __name__ == "__main__":

    # /env/cns/proj/projet_CSD/scratch/assemblies/Ecro_F_AB1/

    parser = argparse.ArgumentParser(description='Run the metagenomic pipeline')
    parser.add_argument('-cwd', help='the path to the working directory where the assembly folder is present')
    parser.add_argument('-collection', help='the name of the project that will be prefixed in the contig names')
    parser.add_argument('-cpu',type=int,default=1,help='number of CPUs used by hhblits (default: 1)')
    args = parser.parse_args()

    # checking arguments

    if args.collection == None :
        sys.exit('the collection parameter is empty, please provide one collection')
    else:
        collection = args.collection

    if not os.path.exists(args.cwd) :
        sys.exit(args.cwd+' does not exist')
    else:
        directory = os.path.abspath(args.cwd)

    if args.cpu < 1 :
        cpu = str(1)
    else:
        cpu = str(args.cpu)


    config_filename = directory+'/'+'assembly'+'/'+'info.json'

    # updating the assembly
    if not os.path.exists(config_filename) :
        update(directory)

    # loading the json info file
    if not os.path.exists(config_filename) :
        sys.exit(config_filename+' does not exist, exit')
    else:
        config_filename = os.path.abspath(config_filename)
        
    with open(config_filename) as f:
        data = json.load(f)

    project = data['project']
    sample = data['sample']
    profileDb_filename = data['anvio_profileDb_filename']
    contigDb_filename = data['anvio_contigDb_filename']
    contig_filename = data['assembly_contig_filename']
    protein_filename = data['assembly_protein_filename']
    anvio_protein_filename = data['anvio_protein_filename']
    bam_filename = data['assembly_bam_filename']
    bai_filename = data['assembly_bai_filename']
    coverage_contigs_filename = data['anvio_coverage_contigs_filename']
    basic_info_contigs_filename = data['anvio_basic_info_contigs_filename']
    gene_taxo_anvio_filename = data['anvio_gene_taxo_anvio_filename']
    taxo_anvio_filename = data['anvio_taxo_filename']


    print('\n')
    print('##############')
    print('# parameters #')
    print('##############')
    print()
    print('Project name: '+project)
    print('Sample name: '+sample)
    print('Collection name: '+collection)
    print('Working directory: '+directory)
    print('Number of CPUs: '+str(cpu))


    print('\n')
    print('############')
    print('# pipeline #')
    print('############')
    print()


    refiningBins_directory = directory+'/'+'refinedBins'
    if os.path.exists(refiningBins_directory) :
        sys.exit(refiningBins_directory+' already exist, please remove it first')
    os.mkdir(refiningBins_directory)


    #################
    # Step 1: ANVIO #
    #################

    # get the contig name from the assembly 
    asmContig2len = dict()
    for record in SeqIO.parse(contig_filename,'fasta') :
        asmContig2len[record.id] = len(record)

    # running anvi-summarize
    anvio_directory = directory+'/'+'refinedBins'+'/'+'ANVIO'
    if os.path.exists(anvio_directory) :
        sys.exit(anvio_directory+' already exist, please remove it first')
    else:
        os.mkdir(anvio_directory)
        anvi_summarize_directory = anvio_directory+'/'+'SAMPLES-SUMMARY'
        cmd = 'source activate anvio-6.2 && anvi-summarize --report-aa-seqs-for-gene-calls -p '+profileDb_filename+' -c '+contigDb_filename+' -o '+anvi_summarize_directory+' -C '+collection+' >/dev/null 2>&1'
        print(cmd)
        status = os.system(cmd)
        print('status: '+str(status))
        if not status == 0:
            sys.exit('something went wrong with anvi-summarize, exit')
        
        # ADD A FILE WITH THE COLLECTIONS NAME
        output = open(anvio_directory+'/'+'collection.txt','w')
        output.write(collection+'\n')
        output.close()
    
    bin_dir = anvio_directory+'/'+'bins'
    if os.path.exists(bin_dir) :
        sys.exit(bin_dir+' already exist, please remove it first')
    os.mkdir(bin_dir)

    # creating a link to the bin contig files
    partialScaffold2bin = defaultdict(set)
    bin2partialScaffold = defaultdict(set)
    scaffold2bin = dict()
    print(bin_dir)
    for root, dirs, files in os.walk(anvi_summarize_directory+'/'+'bin_by_bin', topdown = False):
        for binName in dirs:
            fasta_filename = anvi_summarize_directory+'/'+'bin_by_bin'+'/'+binName+'/'+binName+'-contigs.fa'
            for record in SeqIO.parse(fasta_filename,'fasta') :
                if record.id not in asmContig2len :
                    completeScaffold = record.id.split('_partial_')[0]
                    partialScaffold2bin[completeScaffold].add(binName)
                    bin2partialScaffold[binName].add(record.id)
                    print(record.id+' ('+str(len(record))+' nt) in '+binName+' does not exist in '+contig_filename+'\t'+completeScaffold+' ('+str(asmContig2len[completeScaffold])+' nt)')
                else:
                    scaffold2bin[record.id] = binName
                    if int(len(record)) != int(asmContig2len[record.id]) :
                        print(record.id+' ('+str(len(record))+' nt) ('+str(asmContig2len[record.id])+' nt)')

            if not os.path.exists(bin_dir+'/'+binName+'.fna') and binName not in bin2partialScaffold :
                os.symlink(fasta_filename,bin_dir+'/'+binName+'.fna')
    
    # creating bin files for bins with partial scaffolds
    partialScaffold2seq = dict()
    for record in SeqIO.parse(contig_filename,'fasta') :
        if record.id in partialScaffold2bin :
            partialScaffold2seq[record.id] = record
        else:
            continue

    for binName in bin2partialScaffold :
        print(binName+' contains partial contigs')
        fasta_filename = anvi_summarize_directory+'/'+'bin_by_bin'+'/'+binName+'/'+binName+'-contigs.fa'
        output_filename = bin_dir+'/'+binName+'.fna'
        output = open(output_filename,'w')
        contigSet = set()
        for record in SeqIO.parse(fasta_filename,'fasta') :
            if record.id in bin2partialScaffold[binName] : # if contig is partial
                print(record.id+'\t'+'partial'+str(partialScaffold2bin[completeScaffold]))
                completeScaffold = record.id.split('_partial_')[0]
                if len(partialScaffold2bin[completeScaffold]) == 1 :# if partial scaffold in only one bin
                    print('\tone bin')
                    if completeScaffold not in contigSet :
                        contigSet.add(completeScaffold)
                        SeqIO.write(partialScaffold2seq[completeScaffold],output,'fasta')
                        scaffold2bin[completeScaffold] = binName
                        print('\t\tcopy')
                    else :
                        print('\t\tdont copy')
                        continue
                else : # else go to the unbinned
                    print('\tseveral bins')
                    continue
            else:
                print(record.id+'\t'+'ok')
                SeqIO.write(record,output,'fasta')
        output.close()

    #  creating an unbinned file
    seqList = list()
    for record in SeqIO.parse(contig_filename,'fasta') :
        if record.id in scaffold2bin :
            continue
        if len(record) < 1000 :
            continue
        seqList.append(record)
        scaffold2bin[record.id] = 'Unbinned'

    unbinned_filename = bin_dir+'/'+'Unbinned.fna'
    SeqIO.write(seqList,unbinned_filename,'fasta')

    # parsing the anvio results
    anvio_scaffold2taxonomy = detectingContigTaxonomy(gene_taxo_anvio_filename , taxo_anvio_filename , anvio_protein_filename )
    anvio_scaffold2info = gettingContigInfo(basic_info_contigs_filename, coverage_contigs_filename )


    prokBinSet = set()
    for root, dirs, files in os.walk(bin_dir) :
        for filename in files :
            if filename == 'Unbinned.fna' :
                continue
            if re.search(r'^Euk',filename) :
                continue
            prokBinSet.add(filename)

    if len(prokBinSet) == 0 :
        writingOutputProkaryote( data, refiningBins_directory , anvio_scaffold2info, anvio_scaffold2taxonomy , partialScaffold2bin)
        sys.exit()

    ###################
    # step 2: refineM #
    ###################

    refineM_dir = refiningBins_directory+'/'+'refineM'
    runningRefineM(refiningBins_directory,refineM_dir, bin_dir,contig_filename,bam_filename,cpu)


    ##################
    # step 3: CheckM #
    ##################

    checkm_dir = refiningBins_directory+'/'+'CheckM'
    runningCheckM(checkm_dir,bin_dir,cpu)



    ###################
    # step 4: GTDB-tk #
    ###################

    gtdbtk_dir = refiningBins_directory+'/'+'GTDB-tk'
    runningGTDBtk(gtdbtk_dir,bin_dir,cpu)


    #######################
    # writing the results #
    #######################

    writingOutput( data, refiningBins_directory , anvio_scaffold2info, anvio_scaffold2taxonomy , partialScaffold2bin)

    #######################
    # writing the results #
    #######################

    database = directory+'/'+'refinedBins'+'/'+'REFINEDBINS.db'

    

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
        populate_bins_table(conn , directory+'/'+'refinedBins'+'/'+'output'+'/'+'Bins_summary.tsv' , sample , project)
        print('done')
        # create scaffolds table
        print('create scaffolds table...')
        create_table(conn, sql_create_scaffolds_table)
        populate_scaffolds_table(conn , directory+'/'+'refinedBins'+'/'+'output')
        print('done')
    else:
        print("Error! cannot create the database connection.")


    right = 0o0666
    os.chmod(database, right )
