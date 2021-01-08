#! /usr/bin/python3

import os,sys,re
from collections import defaultdict
from Bio import SeqIO

def runningCheckM(checkm_dir,gene_dir) :
    print('checkM...')
    os.mkdir(checkm_dir)
    os.mkdir(checkm_dir+'/'+'proteins')
    os.mkdir(checkm_dir+'/'+'output')
    
    for root, dirs, files in os.walk(gene_dir, topdown = False):
        for filename in files :
            if re.search(r'_genes.faa',filename) :
                print(filename)
                binName = filename.replace('_genes.faa','')
                if binName == 'Unbinned' :
                    continue
                if re.match(r'Euk',binName) :
                    continue
                print(root+'/'+filename)
                os.symlink(root+'/'+filename,checkm_dir+'/'+'proteins'+'/'+binName+'.faa')
    cmd = 'conda activate checkM-1.1.3 && checkm lineage_wf --genes -t 1 -x faa '+checkm_dir+'/'+'proteins'+' '+checkm_dir+'/'+'output > '+checkm_dir+'/'+'checkm.log 2>&1'
    print(cmd)
    status = os.system(cmd)
    print('status: '+str(status)+'\n')
    if not status == 0 :
        sys.exit('something went wrong with checkm lineage_wf, exit.')


def writtingOutput(genomicOutliers_filename, taxoOutlier_filename, taxoProfile_dir, refineM_scaffold2info, anvio_scaffold2taxonomy, scaffold2bin) :
    
    ###########
    # refineM #
    ###########
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
    file = open(taxoOutlier_filename,'r')
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

    #########
    # ANVIO #
    #########
    print('writting output...')
    output = open('refineM.outpout','w')     
    output.write('scaffold'+'\t'+'bin'+'\t'+'refineM_outlier'+'\t'+'anvio_length'+'\t'+'anvio_gc'+'\t'+'anvio_nb_splits'+'\t'+'anvio_coverage'+'\t'+'anvio_taxonomy'+'\t'+header+'\n')
    for scaffold,binName in scaffold2bin.items() :
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
        scaffold,coverage = line.split('\t')
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


collection = 'raphael_1_20201108'
directory = '/env/cns/proj/projet_CSD/scratch/assemblies/Ecro_F_AB1'

profileDb_filename = directory+'/'+'assembly'+'/'+'anvio'+'/'+'PROFILE.db'
contigDb_filename =  directory+'/'+'assembly'+'/'+'contigs.db'
scaffold_filename =  directory+'/'+'assembly'+'/'+'megahit.contigs.renamed.fa'
protein_filename =   directory+'/'+'assembly'+'/'+'proteins.anvio.tab'
bam_filename = directory+'/'+'assembly'+'/'+'bt2'+'/'+'megahit.contigs.renamed.fa.bam'
bai_filename = directory+'/'+'assembly'+'/'+'bt2'+'/'+'megahit.contigs.renamed.fa.bam.bai'
cpu = str(36)


# anvio_directory = directory+'/'+'refinedBins'+'/'+'ANVIO'
# bin_dir = anvio_directory+'/'+'bins'
# scaffold2bin = dict()
# print(bin_dir)
# for root, dirs, files in os.walk(bin_dir, topdown = False):
#     for filename in files :
#         binName = filename.replace('.fna','')
#         filename = root+'/'+filename
#         for record in SeqIO.parse(filename,'fasta') :
#             scaffold2bin[record.id] = binName
    

# refiningBins_directory = directory+'/'+'refinedBins'
# refineM_dir = refiningBins_directory+'/'+'refineM'
# genomic_dir = refineM_dir+'/'+'genomicProperties'
# stat_dir = genomic_dir+'/'+'stats'
# genomicOutliers_dir = genomic_dir+'/'+'outliers'
# genomicOutliers_filename = genomicOutliers_dir+'/'+'outliers.tsv'
# taxo_dir = refineM_dir+'/'+'taxonomy'
# taxoProfile_dir = taxo_dir+'/'+'profiles'+'/'+'bin_reports'
# taxoOutlier_filename = taxo_dir+'/'+'outliers'+'/'+'taxon_filter.tsv'
# anvio_directory = directory+'/'+'refinedBins'+'/'+'ANVIO'
# gene_dir = anvio_directory+'/'+'genes'
# checkm_dir = refiningBins_directory+'/'+'CheckM'




datatable_dir = directory+'/'+'assembly'+'/'+'datatables'
# if datatables isn't prensent, create it #
if not os.path.exists(datatable_dir) :
    os.mkdir(datatable_dir)


taxo_anvio_filename = datatable_dir+'/'+'taxon_names.txt'
if not os.path.exists(taxo_anvio_filename) :
    cmd = 'source activate anvio-6.2 && anvi-export-table '+contigDb_filename+' --table taxon_names -o '+taxo_anvio_filename
    print(cmd)
    status = os.system(cmd)
    print('status: '+str(status)+'\n')
    if not status == 0 :
        sys.exit('something went wrong with anvi-export-table, exit.')


gene_taxo_anvio_filename = datatable_dir+'/'+'genes_taxonomy.txt'
if not os.path.exists(gene_taxo_anvio_filename) :
    cmd = 'source activate anvio-6.2 && anvi-export-table '+contigDb_filename+' --table genes_taxonomy -o '+gene_taxo_anvio_filename
    print(cmd)
    status = os.system(cmd)
    print('status: '+str(status)+'\n')
    if not status == 0 :
        sys.exit('something went wrong with anvi-export-table, exit.')


basic_info_contigs_filename = datatable_dir+'/'+'contigs_basic_info.txt'
if not os.path.exists(basic_info_contigs_filename) :
    cmd = 'source activate anvio-6.2 && anvi-export-table '+contigDb_filename+' --table contigs_basic_info -o '+basic_info_contigs_filename
    print(cmd)
    status = os.system(cmd)
    print('status: '+str(status)+'\n')
    if not status == 0 :
        sys.exit('something went wrong with anvi-export-table, exit.')


coverage_contigs_filename = datatable_dir+'/'+'contigs_coverage_info.txt'
if not os.path.exists(coverage_contigs_filename) :
    cmd = 'source activate anvio-6.2 && anvi-export-splits-and-coverages -p '+profileDb_filename+' -c '+contigDb_filename+' -o '+datatable_dir+' -O '+'tmp'+' --report-contigs'
    print(cmd)
    status = os.system(cmd)
    print('status: '+str(status)+'\n')
    if not status == 0 :
        sys.exit('something went wrong with anvi-export-splits-and-coverages, exit.')
    os.rename(datatable_dir+'/'+'tmp-COVs.txt',coverage_contigs_filename)
    os.remove(datatable_dir+'/'+'tmp-CONTIGS.fa' )



# anvio_scaffold2taxonomy = detectingContigTaxonomy(gene_taxo_anvio_filename , taxo_anvio_filename , protein_filename )
# anvio_scaffold2info = gettingContigInfo(basic_info_contigs_filename, coverage_contigs_filename )

# writtingOutput( genomicOutliers_filename, taxoOutlier_filename, taxoProfile_dir, anvio_scaffold2info, anvio_scaffold2taxonomy, scaffold2bin)


# if the bam_filename isn't sorted, do #
if not os.path.exists(bai_filename) :
    tmp_bam_filename = bam_filename+'.unsorted'
    os.rename(bam_filename,tmp_bam_filename)
    cmd = 'samtools sort -@ '+str(cpu)+' -o '+bam_filename+' -O BAM '+tmp_bam_filename
    print(cmd)
    status = os.system(cmd)
    print('status: '+str(status))

    # creating the index file
    cmd = 'samtools index '+bam_filename
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

scaffold2bin = dict()
print(bin_dir)
for root, dirs, files in os.walk(anvi_summarize_directory+'/'+'bin_by_bin', topdown = False):
    for binName in dirs:
        # if re.match(r'Euk',binName) :
        #     continue
        fasta_filename = anvi_summarize_directory+'/'+'bin_by_bin'+'/'+binName+'/'+binName+'-contigs.fa'
        for record in SeqIO.parse(fasta_filename,'fasta') :
            scaffold2bin[record.id] = binName


        if not os.path.exists(bin_dir+'/'+binName+'.fna') :
            print(binName+'\t'+fasta_filename)
            os.symlink(fasta_filename,bin_dir+'/'+binName+'.fna')


#  creating an unbinned file
seqList = list()
for record in SeqIO.parse(scaffold_filename,'fasta') :
    if record.id in scaffold2bin :
        continue
    if len(record) < 1000 :
        continue
    seqList.append(record)
    scaffold2bin[record.id] = 'Unbinned'
unbinned_filename = bin_dir+'/'+'Unbinned.fna'
SeqIO.write(seqList,unbinned_filename,'fasta')
    


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


##########
# CheckM #
##########

checkm_dir = refiningBins_directory+'/'+'CheckM'
runningCheckM(checkm_dir,gene_dir)


###########
# GTDB-tk #
###########

# gtdbtk_dir = refiningBins_directory+'/'+'GTDB-tk'
# runningGTDBtk(gtdbtk_dir,gene_dir,cpu)



#######################
# parsing the results #
#######################

anvio_scaffold2taxonomy = detectingContigTaxonomy(gene_taxo_anvio_filename , taxo_anvio_filename , protein_filename )
anvio_scaffold2info = gettingContigInfo(basic_info_contigs_filename, coverage_contigs_filename )

writtingOutput( genomicOutliers_filename, taxoOutlier_filename, taxoProfile_dir, anvio_scaffold2info, anvio_scaffold2taxonomy, scaffold2bin)
