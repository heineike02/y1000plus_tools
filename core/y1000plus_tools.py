# -*- coding: utf-8 -*-
import os 
base_dir = ''
data_processing_dir = ''
#These need to be set after importing the module based on file structure 
#set in std_libraries.py
#I could probably do it automatically with relative paths. 
home_dir = ''
print("home directory is " + home_dir)
base_dir = ''
print("y1000plus_tools dir is " + base_dir ) 
y1000plus_dir = ''
print("y1000plus data dir is " + y1000plus_dir)



import pandas as pd
import numpy as np
# import matplotlib.pyplot as plt
import gffutils

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import SeqIO
#from ete3 import Tree
#ete3 is not officially supported on windows, and so must be loaded via pip: 
# pip install -U https://github.com/etetoolkit/ete/archive/qt5.zip
# ref: https://groups.google.com/forum/#!topic/etetoolkit/6NblSBPij4o
#20181031: got this error message: twisted 18.7.0 requires PyHamcrest>=1.9.0, which is not installed.

# import re
# import math
# import scipy.stats as stats
# import scipy.spatial.distance as spd
# from collections import Counter
import subprocess
# print('I am importing io_library')

# import requests
# from lxml import etree    #parses xml output
# from itertools import product
# import pickle

# Will need to fix the protein sequence routine, but for promoter sequences, as long as the outgroups aren't included in the species subset
# no need to exclude S.Cer and C.Alb
# missing_specs = {'saccharomyces_cerevisiae', 'candida_albicans'}  #outgroups not included (e.g. 'arthrobotrys_oligospora', 'aspergillus_nidulans'  

scaffold_name_change =  {'metschnikowia_dekortum':'+', 
                         'metschnikowia_borealis': '+', 
                         'metschnikowia_continentalis': '+', 
                         'metschnikowia_kamakouana': '+', 
                         'metschnikowia_hibisci': '+', 
                         'metschnikowia_santaceciliae': '+', 
                         'metschnikowia_ipomoeae': '+',
                         'metschnikowia_mauinuiana': '+',
                         'metschnikowia_shivogae': '+'}
        
scaffold_name_change_specs = scaffold_name_change.keys()

# # spec_lookup = {'Klac' : 'Kluyveromyces lactis', 'Scer': 'Saccharomyces cerevisiae', 
# #  'Cgla' : 'Candida glabrata' , 'Ncas': 'Naumovozyma castellii', 
# #  'Sbay' : 'Saccharomyces bayanus', 'Smik': 'Saccharomyces mikatae',
# #  'Lwal' : 'Lachancea waltii', 'Spar' : 'Saccharomyces paradoxus', 
# #  'Lklu' : 'Lachancea kluyverii', 'Dhan': 'Debaryomyces hansenii', 
# #  'Calb' : 'Candida albicans', 'Ylip': 'Yarrowia lipolytica', 
# #  'Sjap' : 'Schizosaccharomyces japonicus' , 'Spom' : 'Schizosaccharomyces pombe' }

# print("data processing dir is " + data_processing_dir )


def build_y1000_species_table():
    #Adds index for orthogroup numbers to 343taxa table and saves the file as 
    #y1000_species_table.csv
    
    y1000_species_table_fname = y1000plus_dir + "y1000_species_table.csv"
    
    #iterate through protein ids and extract numbers for each species.  
    y1000_species_fname = y1000plus_dir + "343taxa_speicies-name_clade-name_color-code.txt"
    y1000_species = pd.read_table(y1000_species_fname, index_col=0)
    
    spec_old_subset_copy = set(y1000_species.old_speceis_names)
    y1000_genename_lookup_fname = y1000plus_dir + "orthomcl_output/orthomcl_SeqIDs_index.txt"
    
    species_index = {} 
    with open(y1000_genename_lookup_fname, 'r') as f: 
        for line in f: 
            y1000_id, spec_gene = line.split(': ')
            spec, gene_id = spec_gene.split('@')
            if spec in spec_old_subset_copy:
                species_index[spec] = y1000_id.split('_')[0]
                spec_old_subset_copy = spec_old_subset_copy - {spec}

    
    y1000_species = y1000_species.assign(spec_og_id = [species_index[spec] for spec in y1000_species['old_speceis_names']])
    
    y1000_species.rename(columns={'old_speceis_names':'old_species_names', 'speceis_names_fig2':'species_names_fig2'}, inplace=True)
    
    y1000_species.to_csv(y1000_species_table_fname)
    
    return y1000_species
    



def id_lookup_generate(y1000_species_subset):
    #species_subset is a dataframe that is a subset of the dataframe created by load_y1000_species
    #for the subset of species we care about, make lookup tables from 
    #ortho_MCL number to gene ID to actual gene name from 343taxa_protein_IDs_index

    #Lookup between species identifiers
    spec_from_old_lookup = dict(zip(y1000_species_subset.old_species_names, y1000_species_subset.original_genome_id))

    #set of species we are interested using old_speceis_names (sic)
    spec_old_subset = set(y1000_species_subset.old_species_names)
    #spec_old_subset = set(y1000_species_subset[y1000_species_subset['original_genome_id'].isin(y1000_species_extra_list)]['old_speceis_names'])

    #filename with full genenames
    fname_id_full = y1000plus_dir + "/orthomcl_output/343taxa_protein_IDs_index.txt"

    #filename with seqIDs
    y1000_genename_lookup_fname = y1000plus_dir + "orthomcl_output/orthomcl_SeqIDs_index.txt"

    #directory to store lookup table csvs
    dir_id_lookups = y1000plus_dir + "id_lookups/"

    print('linking gene_full to gene_id')
    ids_full = {spec : {} for spec in y1000_species_subset.original_genome_id}
    # prev_spec = 'initialize'
   
    with open(fname_id_full, 'r') as f: 
        for line in f: 
            line_sp = line.split('\t')
            spec_old, gene_id = line_sp[1].split('@')
            jgi_test = line_sp[0][0:3]
            if spec_old in spec_old_subset: 
                if spec_old in {'Candida_albicans','Saccharomyces_cerevisiae'}: 
                    gene_full = gene_id
                elif jgi_test == 'jgi':
                    gene_full = line_sp[0].split('|')[-1]
                else: 
                    gene_full = line_sp[0].split(' ')[1].split('=')[1]                
    #             if prev_spec != spec_old:    #could use this to speed up by gradually decreasing the set of species to search
    #                 print(spec_old)
    #             prev_spec = spec_old
                spec = spec_from_old_lookup[spec_old]
                ids_full[spec][gene_id] = gene_full


    print('Linking y1000_id to gene_id')
    ids_y1000 = {spec : {} for spec in y1000_species_subset.original_genome_id}
    with open(y1000_genename_lookup_fname, 'r') as f: 
        for line in f: 
            y1000_id, spec_gene = line.split(': ')
            spec_old, gene_id = spec_gene.split('@')
            if spec_old in spec_old_subset:
                spec = spec_from_old_lookup[spec_old]
                #species_index[spec] = y1000_id.split('_')[0]
                ids_y1000[spec][gene_id.strip()] = y1000_id


    #looping through all species and making CSVs - better to loop through ids_full because it doesn't have S.Cer
    print('Saving lookup tables as .csv files in ' + dir_id_lookups)
    for spec, ids_full_dict in ids_full.items():  
        print(spec)
        gene_lookup = pd.concat([pd.Series(ids_full_dict, name='gene_full'), pd.Series(ids_y1000[spec], name='y1000_id')], axis=1, sort=True)
        gene_lookup.index.name='gene_id'
        gene_lookup.to_csv(dir_id_lookups + spec + '.csv')
    
    return


def make_og_genes_lookup(y1000_id_list, y1000_species_subset):
    #For genes in a y1000_id_list, makes a dictionary of y1000_ids to orthogroups, and of orthogroups to genes in the orthogroups
    #
    #input:
    #    y1000_id_list:  list of y1000_ids that we want orthogroups for
    #
    #output: 
    #    goi_og_lookup:  dictionary from y1000_ids to orthogroup labels
    #    og_genes_lookup: dictionary from orthogroup labels to list of genes in orthogroups
    
    orthogroup_fname = "/home/heineike/genomes/y1000plus/orthomcl_output/orthomcl_clusters.txt"

    #goi_og_set = set(ohnologs_goi_og_genes['OG_index_low']) | set(ohnologs_goi_og_genes['OG_index_high'])
    #low_lookup = dict(zip(ohnologs_goi_og_genes['OG_index_low'],list(ohnologs_goi_og_genes.index)))
    #high_lookup = dict(zip(ohnologs_goi_og_genes['OG_index_high'],list(ohnologs_goi_og_genes.index)))
    
    y1000_species_subset_ogs = set(y1000_species_subset['spec_og_id'])
    y1000_id_set = set(y1000_id_list)
    goi_og_lookup = {}
    og_genes_lookup = {}
    #single_paralog_in_og = []
    with open(orthogroup_fname, 'r') as f: 
        for line in f: 
            og_genes = line.split()
            og = og_genes[0].strip(':')
            genes = og_genes[1:]
            genes_specs = {gene_ind for gene_ind in genes if int(gene_ind.split('_')[0]) in y1000_species_subset_ogs}
            
            found_gois = list(y1000_id_set & genes_specs)
            #May be more than one found gene in an orthogroup
            if len(found_gois) >= 1: 
                
                if len(found_gois)>1: 
                    print('more than one goi in same orthogroup :' + og + ' found_genes: ')
                    print(found_gois)
                for found_goi in found_gois: 
                    goi_og_lookup[found_goi] = og
                og_genes_lookup[og] = genes_specs 
    
    return goi_og_lookup, og_genes_lookup

def make_gtf_dbs(y1000_species_subset): 
    #Make GTF databases for all selected species
    #Only need to do this once

    gtf_dir = y1000plus_dir + "0_332yeast_genomes/332_genome_annotations/gtf/"

    #y1000_species_subset_subset = y1000_species_subset.loc[y1000_species_subset.index>286, :]

    for genome_fname_base in y1000_species_subset['original_genome_id']: 
        #Skipping S.Cerevisiae and Candida Albicans because they aren't set up in the same format
        if genome_fname_base != 'saccharomyces_cerevisiae': 
            print(genome_fname_base)
            gtf_fname = gtf_dir + genome_fname_base + '.max.gtf'
            db_fname = gtf_dir + 'gffutils_dbs/' + genome_fname_base + '.db'

            #Make new database using gffutils
            #I would love for the ID to be the gene name, but CDS and start and stop codons 
            #are the only things annotated.  
            #If that were the case I could use id_spec = {'gene':'gene_id'}
            #
            #I could possibly make that kind of file with GFFread to merge per biostars post:
            #https://www.biostars.org/p/224372/
            #gffread -E merged.gtf -o- > merged.gff3
            #
            #Instead I am using a fancy query to get all items with a certain gene_ID
            #
            #cursor = gtf_db.execute('select * from features where attributes like "%{}%"'.format(gene_name))

            gtf_db = gffutils.create_db(gtf_fname, dbfn=db_fname, force=True, keep_order=True, 
                                        merge_strategy='error', sort_attribute_values=True, 
                                        disable_infer_transcripts=True, disable_infer_genes=True)


            print(genome_fname_base + ' complete')

    return


def extract_promoters(L_prom, og, og_genes, y1000_species_subset, fname_string):
    #for a given orthogroup, extract promoters and put them into a .fasta file
    #L_prom: Length defined for promoter
    
    #Table to look up species by number: 
    genome_name_lookup = dict(zip(y1000_species_subset['spec_og_id'],y1000_species_subset['original_genome_id']))

    promoter_fname = y1000plus_dir + os.path.normpath('promoter_sets/' + og + '_' + fname_string + '.fasta')

    with open(promoter_fname,'w') as f: 

        #for a given set of genes, 
        # group by species
        og_genes_specs = {spec_gene.split('_')[0]: [] for spec_gene in og_genes}
        for y1000_id in og_genes: 
            spec_og_id = y1000_id.split('_')[0]
            og_genes_specs[spec_og_id].append(y1000_id)

        #for each species in the group, gather promoters.  Extract promoter, print promoter as a line in a fasta file
        #>spec y1000_id gene_id gene_full chrm strand start end L
        
        for spec_og_id, genes in og_genes_specs.items(): 
            genome_name = genome_name_lookup[int(spec_og_id)]
            #print(genome_name)         
            #if not(genome_name in missing_specs):  
            
            #load gene_id map based on the species
            gene_lookup_spec_fname = y1000plus_dir + "id_lookups/" + genome_name + '.csv'
            gene_lookup_spec = pd.read_csv(gene_lookup_spec_fname, index_col='y1000_id')
            
            if genome_name=='saccharomyces_cerevisiae': #If S. Cerevisiae, use SGD promoter database
                sc_promoters = pd.read_pickle(base_dir + os.path.normpath('data/Scer_promoters/sc_promoters.pkl'))
                for y1000_id in genes: 
                    gene_id = gene_lookup_spec.loc[y1000_id,'gene_id']
                    prom_seq = sc_promoters.loc[gene_id,:].prom_seq
                    if L_prom>len(prom_seq):
                        print('S.Cerevisiae promoter is only ' + str(len(prom_seq)) + ' bases long, but L_prom=' + str(L_prom))
                    prom_seq_Ltrim = prom_seq[(700-min(700,L_prom)):]
                    sc_common_name = sc_promoters.loc[gene_id,:].sc_common_name
                    f.write('>species=' + genome_name + ' y1000_id=' + y1000_id + ' gene_id=' + gene_id + ' gene_full=' + gene_id+'_'+sc_common_name+ ' L=' + str(len(prom_seq_Ltrim)) + '\n')
                    f.write(prom_seq_Ltrim + '\n')  
            elif genome_name=='candida_albicans': #If C alb, use CGD based promoter database
                ca_promoters = pd.read_pickle(base_dir + os.path.normpath('data/Calb_promoters/Calb_promoters.pkl'))
                for y1000_id in genes: 
                    print(y1000_id)
                    gene_id = gene_lookup_spec.loc[y1000_id,'gene_id']
                    prom_seq = ca_promoters.loc[gene_id,:].prom_seq
                    if L_prom>len(prom_seq):
                        print('C.Albicans promoter is only ' + str(len(prom_seq)) + ' bases long, but L_prom=' + str(L_prom))
                    prom_seq_Ltrim = prom_seq[(1000-min(1000,L_prom)):]
                    f.write('>species=' + genome_name + ' y1000_id=' + y1000_id + ' gene_id=' + gene_id + 'L=' + str(len(prom_seq_Ltrim)) + '\n')
                    f.write(prom_seq_Ltrim + '\n') 
            else: 

                # Load GTF for given sequence: 
                gtf_dir = y1000plus_dir + "0_332yeast_genomes/332_genome_annotations/gtf/"
                db_fname = gtf_dir + 'gffutils_dbs/' + genome_name + '.db'

                gtf_db = gffutils.FeatureDB(db_fname)

                #Extract sequences from genome
                genome_dir = y1000plus_dir + "0_332yeast_genomes/332_genome_assemblies/"

                genome_fname = genome_dir + genome_name + '.fas'

                #For each gene, extract all related features, and then output chromosome and coordinates for the promoter
                for y1000_id in genes: 
                    gene_full = gene_lookup_spec.loc[y1000_id,'gene_full']
                    gene_id = gene_lookup_spec.loc[y1000_id,'gene_id']

                    cursor = gtf_db.execute('select * from features where attributes like "%' + gene_full + '%"')
                    all_features = cursor.fetchall()
                    if len(all_features) == 0:
                        print('No features found ' + gene_full + ' ' +genome_name)

                    strand = all_features[0]['strand']
                    chrom = all_features[0]['seqid']

                    starts = []
                    ends = []
                    for feature in all_features: 
                        starts.append(feature['start'])
                        ends.append(feature['end'])

                    #Adjust coordinates to get L_prom "promoter" sequences
                    if strand == '-': 
                        prom_end = max(ends) + 1
                        prom_start = prom_end + L_prom   #should do min of this and the total length of the scaffold, 
                    elif strand == '+': 
                        prom_end = min(starts) - 1
                        prom_start = max(0,prom_end - L_prom)

                    #extract sequences from genome
                    #Need to reload the iterator each time
                    seq_records = SeqIO.parse(genome_fname, "fasta")
                    scaffold_assigned = False

                    if genome_name in scaffold_name_change_specs:
                        extra_char = scaffold_name_change[genome_name]
                        for seq_record in seq_records:
                            test_scaffold = ''.join(seq_record.id.split(extra_char))
                            if test_scaffold == chrom:
                                scaffold = seq_record
                                scaffold_assigned=True
                    else:     
                        for seq_record in seq_records:
                            test_scaffold = seq_record.id
                            if test_scaffold == chrom:
                                scaffold = seq_record
                                scaffold_assigned=True

                    assert scaffold_assigned, 'Scaffold not assigned for ' + gene_id + ' in ' + genome_name + ', example scaffold: ' + seq_record.id
                    #if strand is negative, check to see if promoter coordinates are at the end of the scaffold

                    #if strand is negative, check to see if promoter coordinates are at the end of the scaffold

                    L_scaffold = len(scaffold)

                    if strand == '-': 
                        if prom_start > L_scaffold: 
                            print('promoter region extends past the scaffold, genome_name = ' + genome_name + ' Gene = ' + gene_id + ', L_prom = ' + str(L_prom))
                            prom_start = L_scaffold
                        if prom_end > L_scaffold: 
                            print('scaffold ends at the end of the gene, genome_name = ' + genome_name + ' Gene = ' + gene_id)
                            prom_end = L_scaffold

                        promoter = scaffold.seq[prom_end:prom_start].reverse_complement()
                    elif strand == '+': 
                        promoter = scaffold.seq[prom_start:prom_end]

                    f.write('>species=' + genome_name + ' y1000_id=' + y1000_id + ' gene_id=' + gene_id + ' gene_full=' + gene_full +
                          ' scaffold=' + chrom + " strand=" + strand + " start=" + str(prom_start) + ' end=' + str(prom_end) +  ' L=' + str(abs(prom_end-prom_start)) + '\n')
                    f.write(str(promoter.upper()) + '\n')  #I wonder why some of the bases were in lower case
    
    return 

def convert_promoters_for_fimo(fname_in):
    #Convert promoter fasta for FIMO, also make dictionary for converting back

    fname_base = fname_in.split('.fasta')[0]
    fasta_out = fname_base + '_fimo.fasta'
    y1000_id_to_gene_id = {}
    y1000_id_to_Lprom = {}
    with open(fasta_out,'w') as f_out: 
        with open(fname_in,'r') as f_in: 
            for line in f_in: 
                line_out = line
                if line[0]=='>':
                    y1000_id = line.split()[1].split('=')[1]
                    line_out = '>' + y1000_id + '\n'
                    species = line.split()[0].split('=')[1]
                    gene_id = line.split()[2].split('=')[1]
                    y1000_id_to_gene_id[y1000_id] = (species, gene_id)
                if line[0]!='>':
                    y1000_id_to_Lprom[y1000_id] = len(line)-1  #subtract one for new line character
                f_out.write(line_out)
    
    return y1000_id_to_gene_id, y1000_id_to_Lprom

def promoter_scan_fimo(promoters_prefix, motif_name, motif_fname, thresh, motif_in_file='All'): 
    
    promoter_fname =  y1000plus_dir + 'promoter_sets/' + promoters_prefix + '_fimo.fasta'
    fname_prefix = promoters_prefix + '_' + motif_name
    output_dir = y1000plus_dir + 'fimo_results' + os.sep

    if motif_in_file == "All":
        motif_arg = []
    else:
        motif_arg = ["--motif",motif_in_file]
    
    fimo_command = ([ home_dir + "meme/bin/fimo",
                      "--oc", output_dir,
                      "--verbosity", "1",
                      "--thresh", str(thresh)] +
                     motif_arg + 
                     [ motif_fname,
                       promoter_fname]
                   )

    print('fimo command:\n' + ' '.join(fimo_command))

    fimo_output = subprocess.run(fimo_command,stdout = subprocess.PIPE) 

    print("fimo output return code = " + str(fimo_output.returncode))

    #change file prefix and delete output other than .txt file
    files_to_change = ['cisml.xml', 'fimo.html', 'fimo.tsv','fimo.xml', 'fimo.gff']

    for file_to_change in files_to_change: 
        full_file_to_change = output_dir +  file_to_change   
        fimo_fname_out = output_dir + fname_prefix + '_' + file_to_change
        os.rename(full_file_to_change, fimo_fname_out)

    fimo_hits = pd.read_table(output_dir + fname_prefix + '_fimo.tsv',
                         engine='python', skipfooter=4)
    
    return fimo_hits

def extract_protein_seqs(og_genes, fname, y1000_species_subset): 
    #Looks up protein sequences for given list of orthogroup genes 
    #
    ## Does not work for S.Cer, C.Alb, outgroup species
    
    proteins_og_fname = y1000plus_dir + os.path.normpath('proteins_og/' + fname + '.fasta')
    
    genome_name_lookup = dict(zip(y1000_species_subset['spec_og_id'],y1000_species_subset['original_genome_id']))
    
    with open(proteins_og_fname,'w') as f: 
        # group by species
        og_genes_specs = {spec_gene.split('_')[0]: [] for spec_gene in og_genes}
        for y1000_id in og_genes: 
            spec_og_id = y1000_id.split('_')[0]
            og_genes_specs[spec_og_id].append(y1000_id)

        for spec_og_id, genes in og_genes_specs.items(): 
            genome_name = genome_name_lookup[int(spec_og_id)]
            print(genome_name)
            #if S.Cer skip finding protein sequences, 

            ## Should add protein sequences in case it is s.cer

            if not(genome_name in missing_specs):    
                #load gene_id map based on the species
                gene_lookup_spec_fname = y1000plus_dir + "id_lookups/" + genome_name + '.csv'
                gene_lookup_spec = pd.read_csv(gene_lookup_spec_fname, index_col='y1000_id')

                #Extract peptide sequences from peptide fasta from genome
                protein_dir = "/home/heineike/genomes/y1000plus/0_332yeast_genomes/332_genome_annotations/pep/"

                protein_fname = protein_dir + genome_name + '.max.pep'

                seq_records = SeqIO.parse(protein_fname, "fasta")

                genes_lookup = gene_lookup_spec.loc[genes]
                genes_lookup_set = set(gene_lookup_spec.loc[genes,'gene_full'].values)

                for seq_record in seq_records:
                    #gene_full = 'augustus_masked-Deha2C-processed-gene-4.36'
                    gene_full = seq_record.description.split()[1].split('=')[1]
                    #print(gene_full)
                    if (gene_full in genes_lookup_set):
                        #find which y1000_id was matched
                        y1000_rlookup = genes_lookup['gene_full'] == gene_full
                        for gene, tf in y1000_rlookup.items(): 
                            if tf:
                                y1000_id=gene
                        gene_id = genes_lookup.loc[y1000_id, 'gene_id']
                        protein_seq = seq_record.seq
                        f.write('>' + genome_name + '_' + gene_id + ' y1000_id=' + y1000_id + ' gene_full=' + gene_full +
                                ' ' + seq_record.description.split()[2] + #this adds in the cds from the original description
                                '\n')
                        f.write(str(protein_seq) + '\n')  #I wonder why some of the bases were in lower case

    return

def seq_key_func(seq): 
    #returns a gene_id string for a given seq object generated from a promoter fasta file
    #one easy option could be gene_id = seq.id, but the id is hidden in the description. 
    description_dict = {item.split('=')[0] : item.split('=')[1] for item in seq.description.split()}
    gene_id = description_dict['species'] + '@' + description_dict['gene_id'] 
    
    return gene_id