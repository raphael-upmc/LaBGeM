#! /usr/bin/python3

import os,sys,re
from Bio import SeqIO
import statistics
from collections import defaultdict
import argparse



# /env/cns/proj/proj_BSI (projet phaeoexplorer)


def renamingContigs(contig_filename,renamed_contig_filename) :
    output = open(renamed_contig_filename,'w')
    file = open(contig_filename,'r')
    for line in file :
        line = line.rstrip()
        if re.match('>',line) :
            header = line.split()[0][1:]
            defline = project+'__'+sample+'__'+header
            output.write('>'+defline+'\n')
        else:
            output.write(line+'\n')
    file.close()
    output.close()



def extractingBam(bam_filename,contig_filename,final_bam_filename,cpu):
    # creating the bed file
    contigSet = set()
    for record in SeqIO.parse(contig_filename,'fasta') :
        contigSet.add(record.id)

    # converting the bam to a sam file
    tmp_sam_filename = bam_filename+'.tmp.sam'
    cmd = 'samtools view -h '+bam_filename+' -o '+tmp_sam_filename
    print(cmd)
    status = os.system(cmd)
    print('status: '+str(status)+'\n')

    tmp_filename = bam_filename+'.tmp.shorter.sam'
    output = open(tmp_filename,'w')
    file = open(tmp_sam_filename,'r')
    for line in file :
        line = line.rstrip()
        liste = line.split('\t')
        if re.match(r'@SQ',line) :
            contig = liste[1].replace('SN:','')
            if contig in contigSet :
                output.write(line+'\n')
            else:
                continue
        elif re.match(r'@',line) :
            output.write(line+'\n')
        else:
            contig = liste[2]
            if contig in contigSet :
                output.write(line+'\n')
            else:
                continue
    output.close()
    file.close()

    # creating and sorting a bam file 
    cmd = 'samtools view -Sb '+' '+tmp_filename+' | samtools sort -@ '+str(cpu)+' -o'+final_bam_filename
    print(cmd)
    status = os.system(cmd)
    print('status: '+str(status)+'\n')

    # creating the index file
    cmd = 'samtools index '+final_bam_filename
    print(cmd)
    status = os.system(cmd)
    print('status: '+str(status)+'\n')

    # removing the two tmp files
    os.remove(tmp_filename)
    os.remove(tmp_sam_filename)


def parsingProdigal(protein_filename,output_filename,contig_filename) :
    contigList = set()
    for record in SeqIO.parse(contig_filename,'fasta') :
        contigList.add(record.id)


    cpt = 0
    output = open(output_filename,'w')
    header = 'gene_callers_id'+'\t'+'contig'+'\t'+'start'+'\t'+'stop'+'\t'+'direction'+'\t'+'partial'+'\t'+'source'+'\t'+'version'+'\t'+'aa_sequence'
    output.write(header+'\n')
    for record in SeqIO.parse(protein_filename,'fasta') :
        liste = record.description.split(' # ')
        orf = str(cpt)
        cpt += 1
        scaffold = '_'.join(liste[0].split('_')[0:-1])
        if scaffold not in contigList :
            continue

        start = str( int(liste[1]) - 1 )
        end = str( int(liste[2]) - 1 )
        if int( liste[-1].split(';')[1].replace('partial=','') ) == 0 :
            partial = '0'
        else:
            partial = '1'

        if liste[3] == '1' :
            strand = 'f'
        else:
            strand = 'r'
        line = orf+'\t'+scaffold+'\t'+start+'\t'+end+'\t'+strand+'\t'+partial+'\t'+'prodigal'+'\t'+'2.6.3'+'\t'+str(record.seq)
        output.write(line+'\n')
    file.close()
    output.close()


def parsingEukRep(euk_filename,prok_filename,contig2split) :
    contig2orig = dict()
    file = open(euk_filename,'r')
    for line in file :
        contig = line.rstrip()
        for split in contig2split[contig] :
            contig2orig[split] = 'EUK'
    file.close()

    file = open(prok_filename,'r')
    for line in file :
        contig = line.rstrip()
        for split in contig2split[contig] :
            contig2orig[split] = 'PROK'
    file.close()
    
    
    for contig,liste in contig2split.items() :
        for split in liste:
            if split not in contig2orig :
                contig2orig[split] = 'UNK'
            else:
                continue

    return contig2orig


def additionalDataTable(annot2data,output_filename,splitList) :

    output = open(output_filename,'w')
    
    categoryList = sorted(list(annot2data.keys()))
    output.write('item_name'+'\t'+'\t'.join( categoryList )+'\n')

    for contig in splitList :
        line = contig
        for category in categoryList :
            result = str(annot2data[category][contig])
            line += '\t'+result
        line += '\n'
        output.write(line)
    output.close()

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Run the metagenomic pipeline')
    parser.add_argument('-fq1', help='the path of the first FASTQ_FILENAME that contains the reads')
    parser.add_argument('-fq2', help='the path of the second FASTQ_FILENAME that contains the reads')
    parser.add_argument('-cwd', help='the path to the working directory where the assembly folder will be created')
    parser.add_argument('-cpu',type=int,default=1,help='number of CPUs used by hhblits (default: 1)')
    parser.add_argument('-k',type=int,default=25000,help='number of contigs to keep for ANVIO (default: 25000)')
    args = parser.parse_args()

    sample = 'test_20201017'
    project = 'phaeoexplorer'

    # checking arguments

    if args.fq1 == None or args.fq2 == None or args.cwd == None :
        sys.exit('fq1, fq2 and cwd are mandatory parameters')

    if not os.path.exists(args.fq1) or not os.path.exists(args.fq2) :
        sys.exit(args.fq1+' or '+args.fq2+' does not exist, exit')
    elif args.fq1 == args.fq2 :
        sys.exit(' fq1 and fq2 are identical, exit')
    else:
        fastq1_filename = os.path.abspath(args.fq1)
        fastq2_filename = os.path.abspath(args.fq2)

    if not os.path.exists(args.cwd) :
        sys.exit(cwd+' does not exist')
    else:
        cwd = os.path.abspath(args.cwd+'/'+'assembly')

    if args.cpu < 1 :
        cpu = 1
    else:
        cpu = args.cpu

    if args.k < 1 :
        k = 25000
    else:
        k = args.k

    print('\n')
    print('##############')
    print('# parameters #')
    print('##############')
    print()
    print('Working directory: '+cwd)
    print('Fastq1: '+fastq1_filename)
    print('Fastq2: '+fastq2_filename)
    print('Number of CPUs: '+str(cpu))
    print('Number of contigs to consider for ANVIO: '+str(k))



    print('\n')
    print('############')
    print('# pipeline #')
    print('############')
    print()


    #################
    # megahit 1.2.9 #
    #################

    contig_filename = cwd+'/'+'contigs.fa'

    print('\n')
    print('Performing the assembly using  megahit...') # megahit will create a directory named assembly
    if not os.path.exists(contig_filename) :
        cmd = 'megahit -1 '+fastq1_filename+' -2 '+fastq2_filename+' -o '+cwd+' --num-cpu-threads '+str(cpu)   ### 1 paired-end library --mem-flag 0
        print(cmd)
        status = os.system(cmd)
        print(status)
        if status == 0 :
            print('creating the working subdirectories '+cwd)
            os.mkdir(cwd+'/'+'annotations')
            os.mkdir(cwd+'/'+'bt2')
            os.mkdir(cwd+'/'+'anvio')
            os.mkdir(cwd+'/'+'taxonomy')
        else:
            sys.exit('something went wrong with megahit, exit')
        print('\n')
    else:
        if not os.path.exists(cwd+'/'+'annotations') :
            os.mkdir(cwd+'/'+'annotations')

        if not os.path.exists(cwd+'/'+'bt2') :
            os.mkdir(cwd+'/'+'bt2')

        if not os.path.exists(cwd+'/'+'anvio') :
            os.mkdir(cwd+'/'+'anvio')

        if not os.path.exists(cwd+'/'+'taxonomy') :
            os.mkdir(cwd+'/'+'taxonomy')

    print('done')





    ########################
    # renaming the contigs #
    ########################

    print('\n')
    print('Renaming the contigs...')
    renamed_contig_filename = cwd+'/'+'contigs.renamed.fa'
    if not os.path.exists(renamed_contig_filename) :
        renamingContigs(contig_filename,renamed_contig_filename) # renaming the contigs + checking for funky characters
    print('done')


    ###########
    # idba-ud #
    ###########

    # print('\n')
    # print('Performing the scaffolding using IDBA-UD...')
    # cmd = 'scaffold '+contig_filename+' '+fastq1_filename+' '+fastq2_filename+' -seed_kmer 100 -min_contig 1000 -o '+cwd
    # print(cmd)
    # print('done')



    #######################
    # mapping with bowtie #
    #######################

    print('\n')
    print('Performing coverage profiling using bowtie2...')
    basename = os.path.basename(renamed_contig_filename)
    bam_filename = cwd+'/'+'bt2'+'/'+basename+'.bam'

    if not os.path.exists(bam_filename) :
        print(bam_filename)
        os.mkdir(output_directory)

        cmd = 'bowtie2-build --threads '+str(cpu)+' '+renamed_contig_filename+' '+cwd+'/'+'bt2'+'/'+basename
        print(cmd)
        status = os.system(cmd)
        print(status)

        cmd = 'bowtie2 -p '+str(cpu)+' -X 1000 -x '+cwd+'/'+'bt2'+'/'+basename+' -1 '+fastq1_filename+' -2 '+fastq2_filename +' | '+'samtools view -Sb >'+bam_filename
        print(cmd)
        status = os.system(cmd)
        print(status)

    print('done')


    #######################################
    # selecting the 25000 longest contigs #
    #######################################

    print('\n')
    print('Selecting the '+str(k)+' longest contigs...')

    lengthList = list()
    for record in SeqIO.parse(renamed_contig_filename,'fasta') :
        lengthList.append(len(record))

    print('\tTotal number of contigs: '+str(len(lengthList)))
    print('\tMedian length: '+str(statistics.median(lengthList)))

    liste = list()
    cpt = 0
    for length in sorted(lengthList,reverse=True) :
        cpt += 1
        liste.append(length)
        if cpt > k :
            break
    print('\tLength of '+str(k)+'th longuest contig: '+str(length))
    print('\tMedian length of the '+str(k)+' longuest contigs: '+str(statistics.median(liste)))

    contig_filename = cwd+'/'+'contigs.renamed'+'.min'+str(length)+'.fa'
    if not os.path.exists(contig_filename) :
        output = open(contig_filename,'w')
        for record in SeqIO.parse(renamed_contig_filename,'fasta') :
            if len(record) > length :
                SeqIO.write(record,output,'fasta')
            else:
                continue
        output.close()
    print('done')


    ####################################################
    # filtering out, indexing and sorting the filename #
    ####################################################

    final_bam_filename = cwd+'/'+'bt2'+'/'+basename+'.min'+str(length)+'.sorted.bam'
    if not os.path.exists(final_bam_filename) :
        extractingBam(bam_filename,contig_filename,final_bam_filename,cpu)



    ############
    # prodigal #
    ############

    print('\n')
    print('Performing the gene prediction using prodigal...')
    protein_filename = cwd+'/'+'proteins.faa'
    if not os.path.exists(protein_filename) :
        cmd = 'prodigal -i '+renamed_contig_filename+' -a '+protein_filename+' -m -p meta >/dev/null 2>/dev/null'
        print(cmd)
        status = os.system(cmd)
        print(status)

    protein_anvio_filename = cwd+'/'+'proteins.anvio.tab'
    if not os.path.exists(protein_anvio_filename) :
        parsingProdigal(protein_filename,protein_anvio_filename,contig_filename) #only on the k contigs
    print('done')



    ######################
    # adding annotations #
    ######################

    # running eukRep
    print('\n')
    print('Adding annotations...')
    print()
    print('\trunning EukRep...')
    eukrep_prok_filename = cwd+'/'+'annotations'+'/'+'eukrepProk.txt'
    eukrep_euk_filename = cwd+'/'+'annotations'+'/'+'eukrepEuk.txt'
    if not os.path.exists(eukrep_euk_filename) :
        cmd = 'EukRep.py -i '+renamed_contig_filename+' -o '+eukrep_euk_filename+' --prokarya '+eukrep_prok_filename+' --seq_names' 
        print(cmd)
        status = os.system(cmd)
        print('status :'+str(status))
        if status != 0 :
            sys.exit()
    print('\tdone')

    print('done')




    #########################
    #########################
    #### ANVIO CMD LINES ####
    #########################
    #########################





    #######################################
    # Creating an anvi’o contigs database # http://merenlab.org/2016/06/22/anvio-tutorial-v2/#creating-an-anvio-contigs-database
    #######################################

    print('\n\nCreating the contig.db file')
    contig_db_filename = cwd+'/'+'contigs.db'
    if not os.path.exists(contig_db_filename) :
        cmd = 'anvi-gen-contigs-database -f '+contig_filename+' -o '+contig_db_filename+' -n '+'\'An example contigs database\''+' --external-gene-calls '+protein_anvio_filename
        print(cmd)
        status = os.system(cmd)
        print(status)

        cmd = 'anvi-run-hmms -c '+contig_db_filename
        print(cmd)
        status = os.system(cmd)
        print(status)
    print('done')


    # anvi-[import|export|show|delete]-misc-data # import annotations
    print('\n')
    print('Importing the annotations into ANVIO...')

    if not os.path.exists(cwd+'/'+'genes_in_splits.txt') :
        cmd = 'anvi-export-table '+contig_db_filename+' --table genes_in_splits -o '+cwd+'/'+'genes_in_splits.txt'
        print(cmd)
        status = os.system(cmd)
        print('status: '+str(status))

    splitList = set()
    contig2split = defaultdict(set)
    file = open(cwd+'/'+'genes_in_splits.txt','r')
    header = next(file)
    for line in file :
        line = line.rstrip()
        liste = line.split('\t')
        split = liste[1]
        contig = '_'.join( split.split('_')[:-2] )
        splitList.add(split)
        contig2split[contig].add(split)
    file.close()

    contig2eukrep = parsingEukRep(eukrep_euk_filename,eukrep_prok_filename,contig2split)

    items_filename = cwd+'/'+'annotations'+'/'+'items_additional_data.txt'
    annot2data = {'EukRep' : contig2eukrep}

    additionalDataTable(annot2data,items_filename,splitList)
    print('done')


    sys.exit()


    ########################
    # anvi-import-taxonomy #
    ########################

    # running kaiju
    print()
    print('Running Kaiju...')
    gene_call_filename = cwd+'/'+'gene-calls.fa'
    kaiju_filename = cwd+'/'+'taxonomy'+'/'+'kaiju.output'
    kaijuTaxon_filename = cwd+'/'+'taxonomy'+'/'+'kaiju-addTaxonNames.output'

    if not os.path.exists(kaijuTaxon_filename) :

        cmd = 'anvi-get-sequences-for-gene-calls -c '+contig_db_filename+' -o '+gene_call_filename
        print(cmd)
        status = os.system(cmd)
        print('status :'+str(status))

        cmd = 'kaiju -t /env/ig/biobank/by-soft/kaiju/1.7.3/i20200525/nodes.dmp -f /env/ig/biobank/by-soft/kaiju/1.7.3/i20200525/kaiju_db_nr_euk.fmi -i '+gene_call_filename+' -o '+kaiju_filename+' -z '+str(cpu)+' -v'
        print(cmd)
        status = os.system(cmd)
        print('status :'+str(status))

        cmd = 'kaiju-addTaxonNames -t /env/ig/biobank/by-soft/kaiju/1.7.3/i20200525/nodes.dmp -n /env/ig/biobank/by-soft/kaiju/1.7.3/i20200525/names.dmp -i '+kaiju_filename+' -o '+kaijuTaxon_filename+' -r superkingdom,phylum,order,class,family,genus,species'+' -v'
        print(cmd)
        status = os.system(cmd)
        print('status :'+str(status))
    print('done')

    cmd = 'anvi-import-taxonomy-for-genes -i '+kaijuTaxon_filename+' -c '+contig_db_filename+' -p kaiju --just-do-it'
    print(cmd)


    ##########################
    # Creating the profile DB #
    ##########################
    profile_filename = cwd+'/'+'anvio'+'/'+'PROFILE.db'
    if not os.path.exists(profile_filename) :
        cmd = 'anvi-profile -i '+final_bam_filename+' -c '+contig_db_filename+' --sample-name \''+project+'__'+sample+'\' --output-dir '+cwd+'/'+'anvio'+' --cluster-contigs --overwrite-output-destinations'
        print(cmd)
        status = os.system(cmd)
        print('status: '+str(status))

        print('\n\n')
        cmd = 'anvi-import-misc-data '+items_filename+' -p '+profile_filename+' --target-data-table items'
        print(cmd)
        status = os.system(cmd)
        print('status: '+str(status))

    #####################
    # anvio interactive #
    #####################

    cmd = 'anvi-interactive --server-only -p '+profile_filename+' -c '+contig_db_filename
    print(cmd)
    sys.exit()


