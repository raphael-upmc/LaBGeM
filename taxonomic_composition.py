#! /usr/bin/env python


import os,sys
from Bio import SeqIO
import argparse


def running_diamond(path_db,path_query,path_out):
    
    print('diamond...')
    
    cmd = 'source /env/cns/proj/agc/scratch/conda/miniconda3/bin/activate anvio-7 && diamond blastp --query '+path_query+'  --db '+path_db+'  --out '+path_out+' --max-target-seqs 20'
    
    #launch diamond
    status = os.system(cmd)
    print(status)
    if not status == 0:
        sys.exit('something went wrong with diamond, exit')
        
            
def max_hits(result_diamond):
    
    gene_id2max = dict()  
    gene_id2percent_identity = dict()
    
    #retrieve the gene id with the maximum score knowing that we have a lot of gene id with the same maximum score.
    file = open(result_diamond,"r")
    for line in file:
        line = line.rstrip()
        liste = line.split('\t')
        gene_id = liste[0]
        percent_identity = float(liste[2])
        score = float(liste[11])       
        
        if gene_id not in gene_id2max:
            gene_id2max[gene_id] = score
            gene_id2percent_identity[gene_id] = percent_identity
            
        else:
            if score > gene_id2max[gene_id]:
                gene_id2max[gene_id] = score
                gene_id2percent_identity[gene_id] = percent_identity    
            else:
                continue
    file.close()
    return gene_id2max,gene_id2percent_identity  
    
            
def max_gene_id2genome_id(result_diamond,gene_id2max):
    
    #Retrieve the gene id unique with list of genomes id if the are the same score
    
    gene_id2gene_db_list = dict()
    genome = list()
    file = open(result_diamond,"r")
    for line in file :
        line = line.rstrip()
        liste = line.split('\t')
        gene_id = liste[0]
        gene_db = liste[1]
        score = float(liste[11])
        maxi = gene_id2max[gene_id]
        if maxi == score:
            if gene_id not in gene_id2gene_db_list :
                gene_id2gene_db_list[gene_id] = []
                gene_id2gene_db_list[gene_id].append(gene_db)
            else:
                gene_id2gene_db_list[gene_id].append(gene_db)
        else:
            continue
    file.close()
    return gene_id2gene_db_list


def gtdb_taxonomy_bank(path_taxonomy):
    
    #create dictionary
    genome_id2taxonomy = dict()
    file = open(path_taxonomy,"r")   
    for line in file:
        line = line.rstrip()
        liste = line.split('\t')
        genome_id = liste[0]
        gtdb_taxonomy = liste[1]
        genome_id2taxonomy[genome_id] = gtdb_taxonomy
            
    file.close()
    return genome_id2taxonomy


def add_taxonomy(gene_id2gene_db_list,genome_id2taxonomy):
    
    #recover the taxonomy of each genomes:
    
    gene_id2taxonomy = dict()
    for gene_id in gene_id2gene_db_list:
        liste = gene_id2gene_db_list[gene_id]
        if len(liste) == 1:
            taxonomy = genome_id2taxonomy[liste[0]]
            gene_id2taxonomy[gene_id] = taxonomy
        else:
            element2taxonomy = dict()
            for genome_id in liste:
                #print(len (liste))
                taxonomy = genome_id2taxonomy[genome_id]
                liste = taxonomy.split(';')
                keys = ['domain','phylum','classe','ordre','famille','genre','espece']
                p = 0
                for key in keys:
                    if key not in element2taxonomy :
                        element2taxonomy[key] = []
                        element2taxonomy[key].append(liste[p])
                        p = p + 1
                    else:
                        element2taxonomy[key].append(liste[p])
                        p = p + 1
                        
            taxonomy2 = list()    
            for key in keys:
                if(len(set(element2taxonomy[key]))==1):
                    taxonomy2.append(element2taxonomy[key][0])            
                else:
                    break
            gene_id2taxonomy[gene_id] = ';'.join(taxonomy2)
    return gene_id2taxonomy
   
        
def gene_id_contigs(result_diamond):
    
    #retrieval gene_id to contig
    
    gene_id2contig_id = dict()
    file = open(result_diamond,"r") 
    for line in file :
        line = line.rstrip()
        liste = line.split('\t')
        gene_id = liste[0]
        gene_id_liste = gene_id.split('_')
        contig_id = '_'.join(gene_id_liste[:-1]) 
        gene_id2contig_id[gene_id] = contig_id
    file.close()
    return gene_id2contig_id


def anvio_coverage(contigs_coverage):
    #retrieval the contigs and Anvio covers
    contig_id2coverage = dict()
    file = open(contigs_coverage,"r")
    next(file)
    for ligne in file :
        ligne = ligne.rstrip()
        liste = ligne.split('\t')
        contig_id = liste[0]
        cover = liste[1]
        contig_id2coverage[ contig_id ] = cover
    file.close()
    return contig_id2coverage


def refineM_coverage(contigs_coverage):
     #retrieval the contigs and refine_M covers
    
    contig_id2coverage_refineM = dict()
    file = open(contigs_coverage,"r")
    for ligne in file :
        ligne = ligne.rstrip()
        liste = ligne.split('\t')
        contig_id = liste[0]
        refineM_cover = liste[2]
        contig_id2coverage_refineM[ contig_id ] = refineM_cover
        
    file.close()
    return contig_id2coverage_refineM


def bin_name(path):
     #retrieval the contigs and bin_name:
    contig2bin = dict()
    filelist = os.listdir(path)
    for bine in filelist:
        bine_name = bine
        with open(path + bine + '/'+ bine+'-contigs.fa', 'r') as fasta_file:
            for record in SeqIO.parse(fasta_file, "fasta"):
                contigs = record.id
                contig2bin[contigs] = bine_name
                
    return contig2bin       
  
       
def kaiju_taxonomy(anvio_file,prodigal_file,kaiju_bank):
    
     #Recuperate anvio_id ---> list(contigs, start ,stop)
    
    anvio_gene_id2liste = dict()
    coord2anvio = dict()
    
    file = open(anvio_file,"r")
    next(file)
    for ligne in file :
        ligne = ligne.rstrip()
        liste = ligne.split('\t')
        anvio_id = liste[0]
        contigs = liste [1]
        start = str(int(liste[2])+1)
        stop = str(int(liste[3])+1)
        anvio_gene_id2liste[anvio_id] = []
        anvio_gene_id2liste[anvio_id].append(contigs)
        anvio_gene_id2liste[anvio_id].append(start)
        anvio_gene_id2liste[anvio_id].append(stop) 
        coord = '_'.join([contigs,start,stop])
        coord2anvio[coord] = anvio_id
    file.close()
    
    #Recuperate prodigale_id ---> list(contigs, start ,stop)
    coord2prodigal = dict()
    prodigal_gene_id2liste = dict()
    file = open(prodigal_file,"r")
    for ligne in file :
        ligne = ligne.rstrip()
        if ligne[0] != '>':
            continue
        liste = ligne.split(' # ')
        prodigal_id = liste[0].replace(">",'')
        start = str(liste[1])
        stop = str(liste[2])
        liste2 = prodigal_id.split('_')
        contigs_id = '_'.join(liste2[:-1]).replace(">",'')
        coord = '_'.join([contigs_id,start,stop])
        coord2prodigal[coord] = prodigal_id
    file.close()
    
    #Recuperate dict prodigale_id --->  anvio_id
    prodigale2anvio = dict()
    for coord , prodigale in coord2prodigal.items():
        if coord in coord2anvio:
            anvio = coord2anvio[coord]
            prodigale2anvio[prodigale] = anvio
        else:
            continue

    #Creation dict anvio_id ----> kaiju_taxonomy
    anvio_id2taxonomy = dict()
    
    file = open(kaiju_bank,"r")
    for ligne in file :
        ligne = ligne.rstrip()
        if ligne[0] != 'C':
            continue
        liste = ligne.split('\t')
        anvio_id = liste [1]
        kaiju_taxonomy = liste[7]
        anvio_id2taxonomy[anvio_id] = kaiju_taxonomy
    file.close()
    return prodigale2anvio,anvio_id2taxonomy


def output_file(gene_id2taxonomy,gene_id2percent_identity,gene_id2max,gene_id2contig_id,contig_id2coverage,contig_id2coverage_refineM,prodigale2anvio,anvio_id2taxonomy,contig2bin):
    
    #writting the output file    
    output = open(path_out+'/'+'result_final_'+marker+'.txt', 'w')
    file = open(result_diamond,"r")
    header = "gene_id	Marker_name	percent_identity	bitscore	anvio_coverage	refineM_coverage	gtdb_taxonomy	kaiju_taxonomy	bin"
    output.write(header+'\n')
    print(header)
    
    for gene_id, taxonomy in gene_id2taxonomy.items():
        marker_name = marker
        percent_identity = gene_id2percent_identity[gene_id]
        bitscore = gene_id2max[gene_id]
        contig_id = gene_id2contig_id[gene_id]
        refineM_coverage = contig_id2coverage_refineM [contig_id]
        #add anvio taxonomy
        if gene_id not in prodigale2anvio:
            anvio_id = 'NA'
        else:
            anvio_id = prodigale2anvio[gene_id]
        #print(anvio_id)
        if anvio_id != 'NA':
            if anvio_id not in anvio_id2taxonomy:
                kaiju_taxonomy = "UNKNOWN"
            else:
                kaiju_taxonomy = anvio_id2taxonomy[anvio_id]
        else: 
            kaiju_taxonomy = "NA"

        if contig_id not in contig_id2coverage:
            coverage = "NA" 
        else:
            coverage = contig_id2coverage[contig_id]
            
        if contig_id not in contig2bin:
            bin_name = "UNBINNED"
        else:
            bin_name = contig2bin[contig_id]
        output.write(gene_id+'\t'+marker_name+'\t'+str(percent_identity)+'\t'+str(bitscore)+'\t'+coverage+'\t'+refineM_coverage+'\t'+taxonomy+'\t'+kaiju_taxonomy+'\t'+bin_name+'\n')
        #print(gene_id+'\t'+marker_name+'\t'+str(percent_identity)+'\t'+str(bitscore)+'\t'+coverage+'\t'+refineM_coverage+'\t'+taxonomy+'\t'+kaiju_taxonomy+'\t'+bin_name+'\n')
    
        
        
if __name__ == "__main__":
    
    
    parser = argparse.ArgumentParser(description='Run the taxonomic composition')
    parser.add_argument('sample', help='the name of the sample that will be analysed')
    parser.add_argument('--markerlist', help='the name of the marker')
    parser.add_argument('--cwd', help='the path to the working directory')   
    args = parser.parse_args()
    
    
    # checking arguments
    
    if args.sample == None :
        sys.exit('please entrer a sample name')
    else:
        sample = args.sample
    
    
    ribosomal_list = ['Ribosomal_L1','Ribosomal_L2','Ribosomal_L3','Ribosomal_L4','Ribosomal_L6','Ribosomal_L9_C','Ribosomal_S11','Ribosomal_S20p','Ribosomal_S2','Ribosomal_S3_C','Ribosomal_S6','Ribosomal_S7','Ribosomal_S8','Ribosomal_S9','Ribosomal_L13','Ribosomal_L16','Ribosomal_L17','Ribosomal_L20','Ribosomal_L21p','Ribosomal_L22','ribosomal_L24','Ribosomal_L27A']
    if args.markerlist == None :
        markerlist = ribosomal_list                  
    else:
        markerlist = args.markerlist.split(',')
        for marker in markerlist:
            if marker not in ribosomal_list:
                sys.exit(marker+' marker not in ribosomal_list, '+str(ribosomal_list))
     
    
    if args.cwd == None:
        sys.exit('please provide cwd option')
    if not os.path.exists(args.cwd) :
        sys.exit(cwd+' does not exist, for phaeoexplorer project enter this cwd : "/env/cns/proj/projet_CSD/scratch/assemblies/" ')
    else:
        cwd = os.path.abspath(args.cwd+'/'+sample)
        if not os.path.exists(cwd):
            sys.exit(cwd+' does not exist')

    path_out = cwd+"/taxonomic_composition/"
    if os.path.exists(path_out) :
        sys.exit(path_out+' exist, please remove taxonomic_composition ')
            
    
    print('## taxonomic composition pipeline ##')
    
    for ribosomal in markerlist:
        marker = ribosomal
        #sample = 'Dher_U_AS1'
        #cwd = "/env/cns/proj/projet_CSD/scratch/assemblies/"
    
        path_db = "/env/cns/proj/agc/scratch/conda/miniconda3/envs/anvio-7/lib/python3.6/site-packages/anvio/data/misc/SCG_TAXONOMY/GTDB/SCG_SEARCH_DATABASES/"+marker+".dmnd"
        path_query = cwd+"/assembly/proteins.faa"
        result_diamond = cwd+"/taxonomic_composition/"+'result_'+marker+'.txt'
    
    
        ### Etape 1 : checking input files
        if not os.path.exists(path_db) :
            sys.exit(path_db+' does not exist')
        else:
            print('path_db is exist..')
    
    
        if not os.path.exists(path_query) :
            sys.exit(path_query+' does not exist')
        else:
            print('path_query is exist..')

      
        ### Etape 2 : checking output file :
        path_out = cwd+"/taxonomic_composition/"
        if not os.path.exists(path_out) :
            os.mkdir(path_out)
    
        if os.path.exists(result_diamond) :
            sys.exit(result_diamond+' already exist, please remove it first')
    
    
        ### Etape 3 :  running diamond
        running_diamond(path_db,path_query,result_diamond)
        print(result_diamond) 
    
   
        ### Etape 4 : recover id genes with max hits
    
        gene_id2max,gene_id2percent_identity = max_hits(result_diamond)
    

        ### Etape 5 : retreive gene id with liste with liste of genome_id
    
        gene_id2gene_db_list = max_gene_id2genome_id(result_diamond,gene_id2max)
    
    
        ### Etape 6 : retreive gtdb taxonomy
    
        path_taxonomy = "/env/cns/proj/agc/scratch_microscope/Data/GTDBtk/707/taxonomy/gtdb_taxonomy.tsv"
        genome_id2taxonomy = gtdb_taxonomy_bank(path_taxonomy)
    
    
        ### Etape 7 : assigne gtdb taxonomy
    
        gene_id2taxonomy = add_taxonomy(gene_id2gene_db_list,genome_id2taxonomy)
    
    
        ### Etape 8 : add gene_id ---> contigs
        gene_id2contig_id = gene_id_contigs(result_diamond)
    
    
        ### Etape 9: add anvio coverage "contigs ---> covers":
        contigs_coverage = cwd+"/assembly/datatables/contigs_coverage_info.txt"
        contig_id2coverage = anvio_coverage(contigs_coverage)
    
    
        ### Etape 10 : add RefineM_coverage "contigs ---> covers":
        contigs_coverage = cwd+"/refinedBins/refineM/genomicProperties/stats/coverage.tsv"
        contig_id2coverage_refineM = refineM_coverage(contigs_coverage)

    
        ### Etape 11 : retrieval the contigs and bin_name
        path= cwd+"/refinedBins/ANVIO/SAMPLES-SUMMARY/bin_by_bin/"
        contig2bin = bin_name(path)
    
    
        ### Etape 12 : add kaiju taxonomy
        path_new = cwd+"/assembly"
        anvio_file = path_new +'/proteins.anvio.tab'
        prodigal_file = path_new +'/proteins.faa'
        kaiju_bank = path_new +'/taxonomy/kaiju-addTaxonNames.output'
    
        prodigale2anvio,anvio_id2taxonomy = kaiju_taxonomy(anvio_file,prodigal_file,kaiju_bank)
    
    
        ### Etape 13 : writting the output file
        result_diamond = cwd+"/taxonomic_composition/"+'result_'+marker+'.txt'
        output_file(gene_id2taxonomy,gene_id2percent_identity,gene_id2max,gene_id2contig_id,contig_id2coverage,contig_id2coverage_refineM,prodigale2anvio,anvio_id2taxonomy,contig2bin)  
  
