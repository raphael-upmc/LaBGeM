#! /usr/bin/env python

import os,sys
from Bio import SeqIO



def diamond_blast(sample,marker):
    #sample = 'Scan_U_DG1'
    path_query = "/env/cns/proj/projet_CSD/scratch/assemblies/"+sample+"/assembly/proteins.faa"
    path_db = "/env/cns/proj/agc/scratch/conda/miniconda3/envs/anvio-7/lib/python3.6/site-packages/anvio/data/misc/SCG_TAXONOMY/GTDB/SCG_SEARCH_DATABASES/"+marker+".dmnd"
    path_out = "/env/cns/proj/projet_CSD/scratch/assemblies/"+sample+"/taxonomic_composition/"
    path_taxonomy = "/env/cns/proj/agc/scratch_microscope/Data/GTDBtk/707/taxonomy"
    
    #print(path_db)

    cmd = 'source /env/cns/proj/agc/scratch/conda/miniconda3/bin/activate anvio-7 && diamond blastp --query '+path_query+'  --db '+path_db+'  --out '+path_out+'result_'+sample+'.txt'+' --max-target-seqs 20'

    #cmd = 'source /env/cns/proj/agc/scratch/conda/miniconda3/bin/activate anvio-7 && diamond blastp --query '+path_query+'  --db '+path_db+'  --out '+path_out+'result_'+sample+'.txt'+' --max-target-seqs 20'
    print(cmd)
    #blast de tous les proteins vs une gene marqueurs:

    ##launch diamond 
    status = os.system(cmd)
    print(status)
    if not status == 0:
        sys.exit('something went wrong with diamond, exit')
        
    gene_id2max = dict()  
    gene_id2evalue = dict()
    gene_id2gene_name = dict()
    #add dict of percent_identity
    gene_id2percent_identity = dict()
    
    file = open(path_out+'/'+'result_'+sample+'.txt',"r")
    ###retreive the geneid with the max score
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
    #print(gene_id2max)
    ### retreive the gene_id unique with liste of genome_id if the are the same score
    gene_id2gene_db_list = dict()
    genome = list()
    file = open(path_out+'/'+'result_'+sample+'.txt',"r")
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
    
    #print(gene_id2gene_db_list)
    #for gene_id,liste in gene_id2gene_db_list.items() :
     #   print(gene_id+'\t'+str(gene_id2max[gene_id])+'\t'+str(liste))
      
     
        
    #create dictionary
    genome_id2taxonomy = dict()
    file = open(path_taxonomy+'/'+'gtdb_taxonomy.tsv',"r")   
    for line in file:
        line = line.rstrip()
        liste = line.split('\t')
        genome_id = liste[0]
        gtdb_taxonomy = liste[1]
        genome_id2taxonomy[genome_id] = gtdb_taxonomy
            
    file.close()
    
    ###recovery of the taxonomy of each of the genomes:
    
    gene_id2taxonomy = dict()
    for gene_id in gene_id2gene_db_list:
        liste = gene_id2gene_db_list[gene_id]
        if len(liste) == 1:
            taxonomy = genome_id2taxonomy[liste[0]]
            gene_id2taxonomy[gene_id] = taxonomy
        else:
            element2taxonomy = dict()
            for genome_id in liste:
                print(len (liste))
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
                #return element2taxonomy
                #print(element2taxonomy)
            taxonomy2 = list()    
            for key in keys:
                if(len(set(element2taxonomy[key]))==1):
                    taxonomy2.append(element2taxonomy[key][0])            
                else:
                    break
            gene_id2taxonomy[gene_id] = ';'.join(taxonomy2)
        #print(gene_id+'\t'+str(gene_id2taxonomy[gene_id]))
        #print(genome_taxonomy)
   

    ##retrieval gene_id to contig
    
    gene_id2contig_id = dict()
    file = open(path_out+'/'+'result_'+sample+'.txt',"r")
    for line in file :
        line = line.rstrip()
        liste = line.split('\t')
        gene_id = liste[0]
        gene_id_liste = gene_id.split('_')
        contig_id = '_'.join(gene_id_liste[:-1]) 
        gene_id2contig_id[gene_id] = contig_id
    #print(gene_id2contig_id)
    file.close()

            
    #####retrieval the contigs and covers
    
    contigs_coverage = "/env/cns/proj/projet_CSD/scratch/assemblies/"+sample+"/assembly/datatables/contigs_coverage_info.txt"
    contig_id2coverage = dict()
    file = open(contigs_coverage,"r")
    next(file)
    for ligne in file :
        ligne = ligne.rstrip()
        liste = ligne.split('\t')
        contig_id = liste[0]
        cover = liste[1]
        contig_id2coverage[ contig_id ] = cover
        #print(contig_id +'\t'+contig_id2coverage[contig_id])
    file.close()
    
    ##### retrieval the contigs and refine_M covers
    
    contigs_coverage2 = "/env/cns/proj/projet_CSD/scratch/assemblies/"+sample+"/refinedBins/refineM/genomicProperties/stats/coverage.tsv"
    contig_id2coverage_refineM = dict()
    file = open(contigs_coverage2,"r")
    for ligne in file :
        ligne = ligne.rstrip()
        liste = ligne.split('\t')
        contig_id = liste[0]
        refineM_cover = liste[2]
        contig_id2coverage_refineM[ contig_id ] = refineM_cover
        
    file.close()
    
    ##retrieval the contigs and bin_name:
    
    path= "/env/cns/proj/projet_CSD/scratch/assemblies/"+sample+"/refinedBins/ANVIO/SAMPLES-SUMMARY/bin_by_bin/"

    contig2bin = dict()
    filelist = os.listdir(path)
    #print(filelist)
    for bine in filelist:
        bine_name = bine
        #print(bine_name)
        with open(path + bine + '/'+ bine+'-contigs.fa', 'r') as fasta_file:
            for record in SeqIO.parse(fasta_file, "fasta"):
                contigs = record.id
                contig2bin[contigs] = bine_name
            
    #print(contig2bin)
    ##################################################  kaiju_taxonomy   ########################################
    
    ###### recuperate anvio_id ---> list(contigs, start ,stop)
    path_new = "/env/cns/proj/projet_CSD/scratch/assemblies/Dher_U_AS1/assembly/"
    
    anvio_gene_id2liste = dict()
    coord2anvio = dict()
    
    file = open(path_new + '/proteins.anvio.tab',"r")
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
    #print(anvio_id+'\t'+str(anvio_gene_id2liste[anvio_id]))
    
    ### recuperate prodigale_id ---> list(contigs, start ,stop)
    coord2prodigal = dict()
    prodigal_gene_id2liste = dict()
    file = open(path_new + '/proteins.faa',"r")
    for ligne in file :
        ligne = ligne.rstrip()
        if ligne[0] != '>':
            continue
        liste = ligne.split(' # ')
        prodigal_id = liste[0].replace(">",'')
        start = str(liste[1])
        stop = str(liste[2])
        #print(star)
        liste2 = prodigal_id.split('_')
        contigs_id = '_'.join(liste2[:-1]).replace(">",'')
        coord = '_'.join([contigs_id,start,stop])
        coord2prodigal[coord] = prodigal_id
        ###remove '>' in the beginning of contigs names
        
       # prodigal_gene_id2liste[prodigal_id] = [contigs_id,start,stop]
    file.close()
    
    ##### recuperate dict prodigale_id --->  anvio_id
    prodigale2anvio = dict()
    for coord , prodigale in coord2prodigal.items():
        #print(coord)
        if coord in coord2anvio:
            anvio = coord2anvio[coord]
            prodigale2anvio[prodigale] = anvio
        else:
            continue

    #####creation dict anvio_id ----> kaiju_taxonomy
    anvio_id2taxonomy = dict()
    
    file = open(path_new + '/taxonomy/kaiju-addTaxonNames.output',"r")
    for ligne in file :
        ligne = ligne.rstrip()
        if ligne[0] != 'C':
            continue
        liste = ligne.split('\t')
        anvio_id = liste [1]
        kaiju_taxonomy = liste[7]
        anvio_id2taxonomy[anvio_id] = kaiju_taxonomy
    file.close()
    #print('la taille' +str(len(anvio_id2taxonomy)))
    #################################### kaiju taxonomy en haut ########################
    
    ### writting the output file
        
    output = open(path_out+'/'+'result_final_'+marker+'.txt', 'w')
    file = open(path_out+'/'+'result_'+sample+'.txt',"r")
    header = "gene_id	Marker_name	percent_identity	bitscore	anvio_coverage	refineM_coverage	gtdb_taxonomy	kaiju_taxonomy	bin"
    output.write(header+'\n')
    print(header)
    
    for gene_id, taxonomy in gene_id2taxonomy.items():
        #print(gene_id)
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
            
        if contig_id not in contig_id2coverage:
            coverage = "NA" 
        else:
            coverage = contig_id2coverage[contig_id]
            
        if contig_id not in contig2bin:
            bin_name = "UNBINNED"
        else:
            bin_name = contig2bin[contig_id]
        #output.write(gene_id+'\t'+marker_name+'\t'+str(percent_identity)+'\t'+str(bitscore)+'\t'+coverage+'\t'+refineM_coverage+'\t'+taxonomy+'\t'+kaiju_taxonomy+'\t'+bin_name+'\n')
        print(gene_id+'\t'+marker_name+'\t'+str(percent_identity)+'\t'+str(bitscore)+'\t'+coverage+'\t'+refineM_coverage+'\t'+taxonomy+'\t'+kaiju_taxonomy+'\t'+bin_name+'\n')
        

#ribosomal_name = ['Ribosomal_L1','Ribosomal_L2','Ribosomal_L3','Ribosomal_L4','Ribosomal_L6','Ribosomal_L9_C','Ribosomal_S11','Ribosomal_S20p','Ribosomal_S2','Ribosomal_S3_C','Ribosomal_S6','Ribosomal_S7','Ribosomal_S8','Ribosomal_S9','Ribosomal_L13','Ribosomal_L16','Ribosomal_L17','Ribosomal_L20','Ribosomal_L21p','Ribosomal_L22','ribosomal_L24','Ribosomal_L27A']      

    path_out = "/env/cns/proj/projet_CSD/scratch/assemblies/"+sample+"/taxonomic_composition/"

    if not os.path.exists(path_out) :
        os.mkdir(path_out)

sample = 'Dher_U_AS1'
marker = 'Ribosomal_L2'
#for marker in ribosomal_name :   
diamond_blast(sample,marker)
