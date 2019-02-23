# -*- coding: utf-8 -*-
import os 
location_input = input("what computer are you on? a = Ben's laptop, b = gpucluster, c = Ben's desktop, d = other")
location_dict = {'a': "C:\\Users\\BMH_work\\", 'b': "/home/heineike/",
                 'c': "C:\\Users\\Ben\\Documents\\", 'd':'you need to add your location to the location_dict'}
home_dir = location_dict[location_input]
print("home directory is " + home_dir)
base_dir = home_dir + os.path.normpath('github/y1000plus_tools/') + os.sep
print("y1000plus_tools dir is " + base_dir ) 
y1000plus_dir_options = {'b':home_dir + os.path.normpath("genomes/y1000plus") + os.sep, 
                         'c': home_dir + os.path.normpath('github/expression_broad_data/expression_data/promoter_phylogenies/y1000plus') + os.sep
                        }
y1000plus_dir = y1000plus_dir_options[location_input]
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
# import subprocess
# print('I am importing io_library')

# import requests
# from lxml import etree    #parses xml output
# from itertools import product
# import pickle

missing_specs = {'saccharomyces_cerevisiae', 'candida_albicans'}  #outgroups not included (e.g. 'arthrobotrys_oligospora', 'aspergillus_nidulans'

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
            #if S.Cer skip finding promoter, 
            if not(genome_name in missing_specs):    
                #load gene_id map based on the species
                gene_lookup_spec_fname = y1000plus_dir + "id_lookups/" + genome_name + '.csv'
                gene_lookup_spec = pd.read_csv(gene_lookup_spec_fname, index_col='y1000_id')

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