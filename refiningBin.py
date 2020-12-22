#! /usr/bin/python3

import os,sys,re
from collections import defaultdict


collection = 'raphael_1_20201108'
directory = '/env/cns/proj/projet_CSD/scratch/assemblies/Ecro_F_AB1'

profileDb_filename = directory+'/'+'assembly/anvio'+'/'+'PROFILE.db'
contigDb_filename = directory+'/'+'assembly'+'/'+'contigs.db'
scaffold_filename = directory+'/'+'assembly'+'/'+'megahit.contigs.renamed.fa'
bam_filename = directory+'/'+'assembly'+'/'+'bt2'+'/'+'megahit.contigs.renamed.fa.bam'

cpu = str(36)


refiningBins_directory = directory+'/'+'refinedBins'
if os.path.exists(refiningBins_directory) :
    print(refiningBins_directory+' already exist, please remove it first')
    sys.exit(refiningBins_directory+' already exist, please remove it first')
os.mkdir(refiningBins_directory)


#########
# ANVIO #
#########

anvio_directory = directory+'/'+'refinedBins'+'/'+'ANVIO'
if os.path.exists(anvio_directory) :
    sys.exit(anvio_directory+' already exist, please remove it first')
else:
    os.mkdir(anvio_directory)
    anvi_summarize_directory = anvio_directory+'/'+'SAMPLES-SUMMARY'
    cmd = 'source activate anvio-6.2 && anvi-summarize -p '+profileDb_filename+' -c '+contigDb_filename+' -o '+anvi_summarize_directory+' -C '+collection
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

print(bin_dir)
for root, dirs, files in os.walk(anvi_summarize_directory+'/'+'bin_by_bin', topdown = False):
    for binName in dirs:
        if re.match(r'Euk',binName) :
            continue
        fasta_filename = anvi_summarize_directory+'/'+'bin_by_bin'+'/'+binName+'/'+binName+'-contigs.fa'
        if not os.path.exists(bin_dir+'/'+binName+'.fna') :
            print(binName+'\t'+fasta_filename)
            os.symlink(fasta_filename,bin_dir+'/'+binName+'.fna')



###########
# refineM #
###########


# Removing contamination based on taxonomic assignments


refineM_dir = refiningBins_directory+'/'+'refineM'
if os.path.exists(refineM_dir) :
    sys.exit(refineM_dir+' already exist, please remove it first')
os.mkdir(refineM_dir)

genomic_dir = refineM_dir+'/'+'genomicProperties'
if os.path.exists(genomic_dir) :
    sys.exit(genomic_dir+' already exist, please remove it first')
os.mkdir(genomic_dir)


stat_dir = genomic_dir+'/'+'stats'
if os.path.exists(stat_dir) :
    sys.exit(stat_dir+' already exist, please remove it first')
os.mkdir(stat_dir)


cmd = 'source activate refineM-0.1.2 && refinem scaffold_stats -c '+cpu+' '+scaffold_filename+' '+bin_dir+' '+stat_dir+' '+bam_filename
print(cmd)
status = os.system(cmd)
print('status: '+str(status))
if not os.path.exists(stat_dir+'/scaffold_stats.tsv'):
    sys.exit('something went wrong with refinem scaffold_stats, exit')


outliers_dir = genomic_dir+'/'+'outliers'
if os.path.exists(outliers_dir) :
    sys.exit(outliers_dir+' already exist, please remove it first')
os.mkdir(outliers_dir)

cmd = 'source activate refineM-0.1.2 && refinem outliers '+stat_dir+'/scaffold_stats.tsv'+' '+outliers_dir
print(cmd)
status = os.system(cmd)
print('status: '+str(status))
if not status == 0:
    sys.exit('something went wrong with refinem scaffold_stats, exit')


filter_dir = genomic_dir+'/'+'filter'
if os.path.exists(filter_dir) :
    sys.exit(filter_dir+' already exist, please remove it first')
os.mkdir(filter_dir)

cmd = 'source activate refineM-0.1.2 && refinem filter_bins '+bin_dir+' '+outliers_dir+'/outliers.tsv '+filter_dir
print(cmd)
status = os.system(cmd)
print('status: '+str(status))
if not status == 0:
    sys.exit('something went wrong with refinem filter_bins, exit')
print('\n\n')


####################
# refineM taxonomy #
####################

gene_dir = anvio_directory+'/'+'genes'
if os.path.exists(gene_dir) :
    sys.exit(gene_dir_dir+' already exist, please remove it first')
os.mkdir(gene_dir)

cmd = 'source activate refineM-0.1.2 && refinem call_genes -c '+cpu+' '+bin_dir+' '+gene_dir
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


cmd = 'source activate refineM-0.1.2 && refinem taxon_profile -c '+cpu+' '+gene_dir+' '+stat_dir+'/scaffold_stats.tsv'+' '+reference_db_filename+' '+reference_taxonomy_filename+' '+taxoProfile_dir
print(cmd)
status = os.system(cmd)
print('status: '+str(status))
if not status == 0:
    sys.exit('something went wrong with refinem taxon_profile, exit')



taxoOutlier_dir = taxo_dir+'/'+'outliers'
if os.path.exists(taxoOutlier_dir) :
    sys.exit(taxoOutlier_dir+' already exist, please remove it first')
os.mkdir(taxoOutlier_dir)

taxonFilter_filename = taxoOutlier_dir+'/'+'taxon_filter.tsv'
cmd = 'source activate refineM-0.1.2 && refinem taxon_filter -c '+cpu+' '+taxoProfile_dir+' '+taxonFilter_filename
print(cmd)
status = os.system(cmd)
print('status: '+str(status))
if not status == 0:
    sys.exit('something went wrong with refinem taxon_filter, exit')



taxoFilter_dir = taxo_dir+'/'+'filters'
if os.path.exists(taxoFilter_dir) :
    sys.exit(taxoFilter_dir+' already exist, please remove it first')
os.mkdir(taxoFilter_dir)

cmd = 'source activate refineM-0.1.2 && refinem filter_bins '+bin_dir+' '+taxonFilter_filename+' '+taxoFilter_dir
print(cmd)
status = os.system(cmd)
print('status: '+str(status))
if not status == 0:
    sys.exit('something went wrong with refinem filter_bins, exit')
