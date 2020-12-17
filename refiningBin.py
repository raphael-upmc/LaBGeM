#! /usr/bin/python3

import os,sys,re
from collections import defaultdict


collection = 'raphael_1_20201108'
directory = '/env/cns/proj/projet_CSD/scratch/assemblies/Ecro_F_AB1'

output_directory = directory+'/'+'assembly/anvio'+'/'+'SAMPLES-SUMMARY'


profileDb_filename = directory+'/'+'assembly/anvio'+'/'+'PROFILE.db'
contigDb_filename = directory+'/'+'assembly'+'/'+'contigs.db'
scaffold_filename = directory+'/'+'assembly'+'/'+'megahit.contigs.renamed.fa'

bam_filename = None
for root, dirs, files in os.walk(directory+'/'+'assembly'+'/'+'bt2', topdown = False):
    for filename in files :
        if re.search(r'^megahit.contigs.renamed.fa.min\d+.sorted.bam$',filename) :
            bam_filename = root+'/'+filename
            print(bam_filename)

if bam_filename == None :    
    sys.exit('bam_filename not found')


cpu = str(36)

if os.path.exists(output_directory) :
    print(output_directory+' already exist, please remove it first')
#    sys.exit(output_directory+' already exist, please remove it first')
else:
    cmd = 'source activate anvio-6.2 && anvi-summarize -p '+profileDb_filename+' -c '+contigDb_filename+' -o '+output_directory+' -C '+collection
    print(cmd)
    status = os.system(cmd)
    print('status: '+str(status))
    if not status == 0:
        sys.exit('something went wrong with anvi-summarize, exit')


bin_dir = directory+'/'+'refinedBins'+'/'+'anvio'+'/'+'bins'
if not os.path.exists(bin_dir) :
    print('creating directory: '+bin_dir)
    os.mkdir(directory+'/'+'refinedBins')
    os.mkdir(directory+'/'+'refinedBins'+'/'+'anvio')
    os.mkdir(bin_dir)

print(bin_dir)
for root, dirs, files in os.walk(output_directory+'/'+'bin_by_bin', topdown = False):
    for binName in dirs:
        if re.match(r'Euk',binName) :
            continue
        fasta_filename = output_directory+'/'+'bin_by_bin'+'/'+binName+'/'+binName+'-contigs.fa'
        if not os.path.exists(bin_dir+'/'+binName+'.fna') :
            print(binName+'\t'+fasta_filename)
            os.symlink(fasta_filename,bin_dir+'/'+binName+'.fna')



# Removing contamination based on taxonomic assignments


refineM_dir = directory+'/'+'refinedBins'+'/'+'anvio'+'/'+'refineM'
if not os.path.exists(refineM_dir) :
    os.mkdir(refineM_dir)

cmd = 'source activate refineM-0.1.2 && refinem scaffold_stats -c '+cpu+' '+scaffold_filename+' '+bin_dir+' '+refineM_dir+' '+bam_filename
print(cmd)
status = os.system(cmd)
print('status: '+str(status))
if not status == 0:
    sys.exit('something went wrong with refinem scaffold_stats, exit')


cmd = 'source activate refineM-0.1.2 && refinem outliers '+refineM_dir+'/scaffold_stats.tsv'+' '+refineM_dir
print(cmd)
status = os.system(cmd)
print('status: '+str(status))
if not status == 0:
    sys.exit('something went wrong with refinem scaffold_stats, exit')


cmd = 'source activate refineM-0.1.2 && refinem filter_bins '+bin_dir+' '+refineM_dir+'/outliers.tsv '+refineM_dir
print(cmd)
status = os.system(cmd)
print('status: '+str(status))
if not status == 0:
    sys.exit('something went wrong with refinem scaffold_stats, exit')


print('\n\n')


gene_dir = directory+'/'+'refinedBins'+'/'+'anvio'+'/'+'genes'
if not os.path.exists(gene_dir) :
    os.mkdir(gene_dir)

cmd = 'source activate refineM-0.1.2 && refinem call_genes -c '+cpu+' '+bin_dir+' '+gene_dir
print(cmd)
status = os.system(cmd)
print('status: '+str(status))
if not status == 0:
    sys.exit('something went wrong with refinem call_genes, exit')



reference_taxonomy_filename = '/env/cns/proj/agc/scratch/proj/GTDB/gtdb_r95_taxonomy.2020-07-30.tsv'
reference_db_filename = '/env/cns/proj/agc/scratch/proj/GTDB/gtdb_r95_protein_db.2020-07-30.faa.dmnd'

taxo_dir = directory+'/'+'refinedBins'+'/'+'anvio'+'/'+'taxonomy'
if not os.path.exists(taxo_dir) :
    os.mkdir(taxo_dir)

cmd = 'source activate refineM-0.1.2 && refinem taxon_profile -c '+cpu+' '+gene_dir+' '+refineM_dir+'/scaffold_stats.tsv'+' '+reference_db_filename+' '+reference_taxonomy_filename+' '+taxo_dir
print(cmd)
status = os.system(cmd)
print('status: '+str(status))
if not status == 0:
    sys.exit('something went wrong with refinem taxon_profile, exit')


taxonFilter_filename = taxo_dir+'/'+'taxon_filter.tsv'
cmd = 'source activate refineM-0.1.2 && refinem taxon_filter -c '+cpu+' '+taxo_dir+' '+taxonFilter_filename
print(cmd)
status = os.system(cmd)
print('status: '+str(status))
if not status == 0:
    sys.exit('something went wrong with refinem taxon_filter, exit')

