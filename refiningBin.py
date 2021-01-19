#! /usr/bin/python3

import os,sys,re
from collections import defaultdict
from Bio import SeqIO
import argparse
import ast
import json
from xlsxwriter.workbook import Workbook
import csv


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


    cmd = 'source activate refineM-0.1.2 && refinem taxon_profile -c '+cpu+' '+gene_dir+' '+stat_dir+'/scaffold_stats.tsv'+' '+reference_db_filename+' '+reference_taxonomy_filename+' '+taxoProfile_dir+' >/dev/null 2>&1'
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
    cmd = 'samtools sort -@ '+str(cpu)+' -o '+bam_filename+' -O BAM '+tmp_bam_filename+' >/dev/null 2>&1'
    print(cmd)
    status = os.system(cmd)
    print('status: '+str(status))

    # creating the index file
    cmd = 'samtools index '+bam_filename+' >/dev/null 2>&1'
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
    output.write('\t\"directory\":\"'+directory+'/'+'assembly'+'\",\n')
    output.write('\t\"contig_filename\":\"'+directory+'/'+'assembly'+'/'+'megahit.contigs.renamed.fa'+'\",\n')
    output.write('\t\"protein_filename\":\"'+directory+'/'+'assembly'+'/'+'proteins.anvio.tab'+'\",\n')
    output.write('\t\"bam_filename\":\"'+bam_filename+'\",\n')
    output.write('\t\"bai_filename\":\"'+bai_filename+'\",\n')
    output.write('\t\"contigDb_filename\":\"'+directory+'/'+'assembly'+'/'+'contigs.db'+'\",\n')
    output.write('\t\"profileDb_filename\":\"'+directory+'/'+'assembly'+'/'+'anvio'+'/'+'PROFILE.db'+'\",\n')
    output.write('\t\"coverage_contigs_filename\":\"'+coverage_contigs_filename+'\",\n')
    output.write('\t\"basic_info_contigs_filename\":\"'+basic_info_contigs_filename+'\",\n')
    output.write('\t\"gene_taxo_anvio_filename\":\"'+gene_taxo_anvio_filename+'\",\n')
    output.write('\t\"taxo_anvio_filename\":\"'+taxo_anvio_filename+'\"\n')
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


def writingOutput(json_data, refiningBins_directory , anvio_scaffold2info, anvio_scaffold2taxonomy, bin2scaffold) :

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
        print(filename)
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
    print(bin2gtdb)



    ###############
    # collections #
    ###############

    filenameList = list()
    output_dir = refiningBins_directory+'/'+'output' 
    os.mkdir(output_dir)

    print('writting outputs in'+output_dir+'...')
    filenameList.append(output_dir+'/'+'Collection.tsv')
    filenameList.append(output_dir+'/'+'Anvio_summary.tsv')
    filenameList.append(output_dir+'/'+'CheckM.tsv')
    filenameList.append(output_dir+'/'+'GTDBtk.tsv')


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
    output = open(output_dir+'/'+'Collection.tsv','w')
    output.write('Project: '+'\t'+json_data['project']+'\n')
    output.write('Sample: '+'\t'+json_data['sample']+'\n')
    output.write('Collection: '+'\t'+collection+'\n')
    output.write('Assembly directory: '+'\t'+json_data['directory']+'\n')
    output.write('Working directory: '+'\t'+refiningBins_directory+'\n')
    output.write('\n')

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
            if binName in bin2checkm and header in bin2gtdb[binName] :
                output.write('\t'+str(bin2gtdb[binName][header]))
            else :
                output.write('\t'+'Na')

        output.write('\n')
    output.close()


    ##########################
    # creating the bin files #
    ##########################
    ###########
    # refineM
    genomicOutliers_filename = refiningBins_directory+'/refineM/genomicProperties/outliers/outliers.tsv'
    taxoOutliers_filename = refiningBins_directory+'/refineM/taxonomy/outliers/taxon_filter.tsv'
    taxoProfile_dir = refiningBins_directory+'/refineM/taxonomy/profiles/bin_reports'

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
                    refineM_scaffold2info[ scaffold ] = line
                file.close()


    ############
    # bin files
    for binName,scaffold2length in bin2scaffold.items() :
        output_filename = output_dir+'/'+binName+'.tsv'
        filenameList.append(output_filename)
        output = open(output_filename,'w')
        output.write('scaffold'+'\t'+'bin'+'\t'+'refineM_outlier'+'\t'+'anvio_length'+'\t'+'anvio_gc'+'\t'+'anvio_nb_splits'+'\t'+'anvio_coverage'+'\t'+'anvio_taxonomy'+'\t'+header+'\n')
        for scaffold,length in sorted(scaffold2length.items() , key = lambda x:x[1] , reverse = True ) :
            if scaffold in anvio_scaffold2taxonomy :
                taxonomy = anvio_scaffold2taxonomy[scaffold]
            else:
                taxonomy = 'Na'

            if scaffold in anvio_scaffold2info :
                info = '\t'.join(anvio_scaffold2info[ scaffold ])
            else:
                info = 'Na\tNa\tNa\tNa'

            if scaffold in refineM_scaffold2info :
                refineM_info = refineM_scaffold2info[scaffold]
            else:
                refineM_info = 'Na'
        
            if scaffold in scaffold2outlier :
                outlier = scaffold2outlier[scaffold]
            else :
                outlier = '-'
            output.write(scaffold+'\t'+binName+'\t'+outlier+'\t'+info+'\t'+taxonomy+'\t'+refineM_info+'\n')
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
        sys.exit(cwd+' does not exist')
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
    profileDb_filename = data['profileDb_filename']
    contigDb_filename = data['contigDb_filename']
    contig_filename = data['contig_filename']
    protein_filename = data['protein_filename']
    bam_filename = data['bam_filename']
    bai_filename = data['bai_filename']
    coverage_contigs_filename = data['coverage_contigs_filename']
    basic_info_contigs_filename = data['basic_info_contigs_filename']
    gene_taxo_anvio_filename = data['gene_taxo_anvio_filename']
    taxo_anvio_filename = data['taxo_anvio_filename']



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
        print(refiningBins_directory+' already exist, please remove it first')
        sys.exit(refiningBins_directory+' already exist, please remove it first')
    os.mkdir(refiningBins_directory)


    #################
    # Step 1: ANVIO #
    #################

    # running anvi-summarize
    anvio_directory = directory+'/'+'refinedBins'+'/'+'ANVIO'
    if os.path.exists(anvio_directory) :
        sys.exit(anvio_directory+' already exist, please remove it first')
    else:
        os.mkdir(anvio_directory)
        anvi_summarize_directory = anvio_directory+'/'+'SAMPLES-SUMMARY'
        cmd = 'source activate anvio-6.2 && anvi-summarize -p '+profileDb_filename+' -c '+contigDb_filename+' -o '+anvi_summarize_directory+' -C '+collection+' >/dev/null 2>&1'
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
    bin2scaffold = dict()
    scaffold2bin = dict()
    print(bin_dir)
    for root, dirs, files in os.walk(anvi_summarize_directory+'/'+'bin_by_bin', topdown = False):
        for binName in dirs:
            if binName not in bin2scaffold :
                bin2scaffold[ binName ] = dict()

            fasta_filename = anvi_summarize_directory+'/'+'bin_by_bin'+'/'+binName+'/'+binName+'-contigs.fa'
            for record in SeqIO.parse(fasta_filename,'fasta') :
                scaffold2bin[record.id] = binName
                bin2scaffold[ binName ][ record.id ] = len(record)

            if not os.path.exists(bin_dir+'/'+binName+'.fna') :
                print(binName+'\t'+fasta_filename)
                os.symlink(fasta_filename,bin_dir+'/'+binName+'.fna')

    #  creating an unbinned file
    bin2scaffold[ 'Unbinned' ] = dict()
    seqList = list()
    for record in SeqIO.parse(contig_filename,'fasta') :
        if record.id in scaffold2bin :
            continue
        if len(record) < 1000 :
            continue
        seqList.append(record)
        scaffold2bin[record.id] = 'Unbinned'
        bin2scaffold[ 'Unbinned' ][ record.id ] = len(record)

    unbinned_filename = bin_dir+'/'+'Unbinned.fna'
    SeqIO.write(seqList,unbinned_filename,'fasta')

    # parsing the anvio results
    anvio_scaffold2taxonomy = detectingContigTaxonomy(gene_taxo_anvio_filename , taxo_anvio_filename , protein_filename )
    anvio_scaffold2info = gettingContigInfo(basic_info_contigs_filename, coverage_contigs_filename )




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

    writingOutput( data, refiningBins_directory , anvio_scaffold2info, anvio_scaffold2taxonomy, bin2scaffold)
