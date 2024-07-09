# -*- coding: utf-8 -*-
### means I need to import the library to my environment
import os 
#base_dir = ''
data_processing_dir = ''
#These need to be set after importing the module based on file structure 
#set in std_libraries.py
#base_dir = ''
#print("y1000plus_tools dir is unset") 
y1000plus_dir = ''
print("y1000plus data dir is unset")

import sys
#import yeast_esr_exp    ## Some functions require this library - especially to look up yeast orfs and common names.    Should break that funciton out separately.  

#These directories are unset at first - they get set externally
# yeast_esr_exp.base_dir = yeast_esr_exp_path 
# yeast_esr_exp.data_processing_dir = yeast_esr_exp_path + os.path.normpath('expression_data') + os.sep



import copy
import pandas as pd
import numpy as np

import matplotlib.cm as cm
import matplotlib.colors as colors
# import matplotlib.pyplot as plt

from statsmodels.distributions.empirical_distribution import ECDF

#import gffutils   #Some functions require this library - it is not in Conda's package list so need to add to environment separately

from Bio.Seq import Seq
#from Bio.Alphabet import generic_dna, IUPAC
from Bio import SeqIO
from Bio import motifs

from ete3 import Tree, SeqMotifFace, TreeStyle, add_face_to_node, RectFace, NodeStyle, TextFace, AttrFace
#from ete3 import Tree
#ete3 is not officially supported on windows, and so must be loaded via pip: 
# pip install -U https://github.com/etetoolkit/ete/archive/qt5.zip
# ref: https://groups.google.com/forum/#!topic/etetoolkit/6NblSBPij4o
#20181031: got this error message: twisted 18.7.0 requires PyHamcrest>=1.9.0, which is not installed.

# import re
# import math
# import scipy.stats as stats
# import scipy.spatial.distance as spd
from collections import OrderedDict  #, Counter
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
    
    y1000_species_table_fname = y1000plus_dir + os.path.normpath("y1000_plus_tools/y1000_plus_tools_data/y1000_species_table.csv")
    
    #iterate through protein ids and extract numbers for each species.  
    y1000_species_fname = y1000plus_dir + 'shen_2018_data' + os.sep + "343taxa_speicies-name_clade-name_color-code.txt"
    y1000_species = pd.read_table(y1000_species_fname, index_col=0)
    
    spec_old_subset_copy = set(y1000_species.old_speceis_names)
    y1000_genename_lookup_fname = y1000plus_dir + os.path.normpath("orthomcl_output/orthomcl_SeqIDs_index.txt")
    
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

def y1000_clade_color_lookup(y1000_species_fname = y1000plus_dir + "y1000_species_table.csv"): 
    y1000_species = pd.read_csv(y1000_species_fname)
    
    clade_color_lookup = {}

    for clade, color in zip(y1000_species['Major clade'], y1000_species['hex']): 
        clade_color_lookup[clade] = color 
        
    return clade_color_lookup

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
    
    orthogroup_fname = y1000plus_dir + os.path.normpath("shen_2018_data/orthomcl_output/orthomcl_clusters.txt")

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

def combine_ogs(og_list, og_genes_lookup): 
    #takes list of ogs and og_genes_lookup and outputs a combined orthogroup name
    #and combined orthogroup list
    ogcomb_name = '_'.join(og_list)

    ogcomb_genes = []
    for og in og_list: 
        ogcomb_genes = ogcomb_genes + list(og_genes_lookup[og])

    return ogcomb_name, ogcomb_genes

def make_gtf_dbs(y1000_species_list): 
    #Make GTF databases for all selected species
    #Only need to do this once
    
    ## Changed input to be list of species names.  Used to be y1000_species_subset
    #will need the 'original_genome_id' column, e.g. y1000_species_subset['original_genome_id']   
    
    gtf_dir = y1000plus_dir + os.path.normpath("shen_2018_data/0_332yeast_genomes/332_genome_annotations/gtf") + os.sep
    
    db_dir = y1000plus_dir + os.path.normpath("y1000plus_tools_data/y1000plus/gffutils_dbs") + os.sep

    #y1000_species_subset_subset = y1000_species_subset.loc[y1000_species_subset.index>286, :]

    for genome_fname_base in y1000_species_list: 
        #Skipping S. cerevisiae and Candida albicans because they aren't set up in the same format
        if genome_fname_base == 'candida_albicans': 
            print(genome_fname_base + ' database not created - format not compatible')
                  
        else:
            print('making ' + genome_fname_base + ' database')
            if genome_fname_base =='saccharomyces_cerevisiae': 
                gtf_fname = y1000plus_dir +  os.path.normpath("y1000plus_tools_data/scer_20181114/saccharomyces_cerevisiae_R64-2-2_20170117.gff") 
                db_fname = y1000plus_dir +  os.path.normpath("y1000plus_tools_data/scer_20181114/saccharomyces_cerevisiae_R64-2-2_20170117.db") 
            else: 
                gtf_fname = gtf_dir + genome_fname_base + '.max.gtf'
                db_fname = db_dir + genome_fname_base + '.db'
            print(gtf_fname)
            print(db_fname)
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

            gtf_db = gffutils.create_db(gtf_fname, 
                                        dbfn=db_fname, 
                                        force=True, 
                                        keep_order=True, 
                                        merge_strategy='error', 
                                        sort_attribute_values=True, 
                                        disable_infer_transcripts=True, 
                                        disable_infer_genes=True)


            print(genome_fname_base + ' complete')

    return


def extract_promoters(L_prom, og, og_genes, y1000_species_subset, fname_string):
    ## Fix error message for short promoters
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
                print('Missing File need to reset file location')
                #sc_promoters = pd.read_pickle(base_dir + os.path.normpath('data/Scer_promoters/sc_promoters.pkl'))
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
                print('Missing File need to reset file location')
                #ca_promoters = pd.read_pickle(base_dir + os.path.normpath('data/Calb_promoters/Calb_promoters.pkl'))
                for y1000_id in genes: 
                    print(y1000_id)
                    gene_id = gene_lookup_spec.loc[y1000_id,'gene_id']
                    prom_seq = ca_promoters.loc[gene_id,:].prom_seq
                    if L_prom>len(prom_seq):
                        print('C.Albicans promoter is only ' + str(len(prom_seq)) + ' bases long, but L_prom=' + str(L_prom))
                    prom_seq_Ltrim = prom_seq[(1000-min(1000,L_prom)):]
                    f.write('>species=' + genome_name + ' y1000_id=' + y1000_id + ' gene_id=' + gene_id + ' L=' + str(len(prom_seq_Ltrim)) + '\n')
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

def extract_all_promoters_for_spec(genome_name, L_prom, fname_out):
    #for a given orthogroup, extract promoters and put them into a .fasta file
    #L_prom: Length defined for promoter

    #Table to look up species by number: 
    #genome_name_lookup = dict(zip(y1000_species_subset['spec_og_id'],y1000_species_subset['original_genome_id']))

    #promoter_fname = y1000plus_dir + os.path.normpath('promoter_sets/' + og + '_' + fname_string + '.fasta')

    with open(fname_out,'w') as f: 

        #load gene_id map based on the species
        gene_lookup_spec_fname = y1000plus_dir + "id_lookups/" + genome_name + '.csv'
        gene_lookup_spec = pd.read_csv(gene_lookup_spec_fname, index_col='y1000_id')

        if genome_name=='saccharomyces_cerevisiae': #If S. Cerevisiae, use SGD promoter database
            print('already have curated promoters for ' + genome_name)
        elif genome_name=='candida_albicans':
            print('already have curated promoters for ' + genome_name)
        else: 
            # Load GTF for given sequence: 
            gtf_dir = y1000plus_dir + "0_332yeast_genomes/332_genome_annotations/gtf/"
            db_fname = gtf_dir + 'gffutils_dbs/' + genome_name + '.db'

            gtf_db = gffutils.FeatureDB(db_fname)

            #Extract sequences from genome
            genome_dir = y1000plus_dir + "0_332yeast_genomes/332_genome_assemblies/"

            genome_fname = genome_dir + genome_name + '.fas'

            Ntot = len(gene_lookup_spec)
            jj = 0

            for y1000_id, (gene_id,gene_full) in gene_lookup_spec.iterrows():
                if np.mod(jj,1000)==0:
                    print(str(jj) + ' of ' + str(Ntot))
                #For each gene, extract all related features, and then output chromosome and coordinates for the promoter

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
                    if prom_end-L_prom<0: 
                        prom_start=0
                        print('Promoter region extendspast the scaffold, genome_name = ' + genome_name + ' Gene = ' + gene_id + ', L_prom = ' + str(L_prom) + ' length avail=' + str(prom_end)) 
                    else: 
                        prom_start = prom_end-L_prom      
                    #prom_start = max(0,prom_end - L_prom)  #taking the max here ensures that you don't go off the end of the scaffold

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

                L_scaffold = len(scaffold)

                if strand == '-': 
                    if prom_start > L_scaffold: 
                        prom_start = L_scaffold
                        print('promoter region extends past the scaffold, genome_name = ' + genome_name + ' Gene = ' + gene_id + ', L_prom = ' + str(L_prom) + ' length avail=' + str(prom_start-prom_end))
                    if prom_end > L_scaffold: 
                        print('scaffold ends at the end of the gene, genome_name = ' + genome_name + ' Gene = ' + gene_id)
                        prom_end = L_scaffold

                    promoter = scaffold.seq[prom_end:prom_start].reverse_complement()
                elif strand == '+': 
                    promoter = scaffold.seq[prom_start:prom_end]

                f.write('>' + genome_name + '@' + gene_id + ' y1000_id=' + y1000_id + ' gene_full=' + gene_full +
                      ' scaffold=' + chrom + " strand=" + strand + " start=" + str(prom_start) + ' end=' + str(prom_end) +  ' L=' + str(len(promoter)) + '\n')
                f.write(str(promoter.upper()) + '\n')  #I wonder why some of the bases were in lower case

                jj=jj+1
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

def promoter_scan_fimo(promoters_fname, fname_prefix, motif_name, motif_fname, thresh, motif_in_file='All'): 
    #promoters_fname: file where promoters are stored
    #fname_prefix: 
    #motif_name: Name of motif that we are searching (will be used to store output)
    #motif_fname: filename of the motif matrix
    #thresh: threshold to use to call a hit. 
    
    meme_dir = '??? need to set'
    
    
    output_dir = y1000plus_dir + 'fimo_results' + os.sep

    if motif_in_file == "All":
        motif_arg = []
    else:
        motif_arg = ["--motif",motif_in_file]
    
    fimo_command = ([ meme_dir + "/bin/fimo",
                      "--oc", output_dir,
                      "--verbosity", "1",
                      "--thresh", str(thresh)] +
                     motif_arg + 
                     [ motif_fname,
                       promoters_fname]
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

def extract_protein_seqs(og_genes, fname_base, y1000_species_subset): 
    #Looks up protein sequences for given list of orthogroup genes 
    #
    ## Does not work for outgroup species
    
    os.mkdir(y1000plus_dir + os.path.normpath('y1000plus_tools_data/y1000plus/proteins_og/' + fname_base))
    proteins_og_fname = y1000plus_dir + os.path.normpath('y1000plus_tools_data/y1000plus/proteins_og/' + fname_base + '/' + fname_base + '.fasta')
    
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
            #if S.Cer or C. Albicans, do slightly different routine
  
            if genome_name == 'saccharomyces_cerevisiae':
                #load gene_id map based on the species
                gene_lookup_spec_fname = y1000plus_dir + os.path.normpath("y1000plus_tools_data/y1000plus/id_lookups/" + genome_name + '.csv')
                gene_lookup_spec = pd.read_csv(gene_lookup_spec_fname, index_col='y1000_id')
                
                protein_fname = y1000plus_dir + os.path.normpath('shen_2018_data/0_332yeast_genomes/332_genome_annotations/Saccharomyces_cerevisiae_S288C_protein.fasta')

                seq_records = SeqIO.parse(protein_fname, "fasta")

                genes_lookup = gene_lookup_spec.loc[genes]
                genes_lookup_set = set(gene_lookup_spec.loc[genes,'gene_id'].values)  #SC specific

                for seq_record in seq_records:
                    gene_id = seq_record.description.split()[0] #SC specific
                    #print(gene_full)
                    if (gene_id in genes_lookup_set):
                        #find which y1000_id was matched
                        y1000_rlookup = (genes_lookup['gene_id'] == gene_id)
                        for gene, tf in y1000_rlookup.items(): 
                            if tf:
                                y1000_id=gene
                        protein_seq = seq_record.seq
                        f.write('>' + genome_name + '_' + gene_id + ' y1000_id=' + y1000_id + ' gene_id=' + gene_id +
                                ' common_name=' + seq_record.description.split()[1] + #this adds in the cds from the original description
                                '\n')
                        f.write(str(protein_seq) + '\n')  #I wonder why some of the bases were in lower case
            elif genome_name == 'candida_albicans':
                #load gene_id map based on the species
                gene_lookup_spec_fname = y1000plus_dir + os.path.normpath("y1000plus_tools_data/y1000plus/id_lookups/" + genome_name + '.csv')
                gene_lookup_spec = pd.read_csv(gene_lookup_spec_fname, index_col='y1000_id')

                protein_fname = y1000plus_dir + os.path.normpath('shen_2018_data/0_332yeast_genomes/332_genome_annotations/Candida_albicans_SC5314_A22_current_default_protein.fasta')
                seq_records = SeqIO.parse(protein_fname, "fasta")

                genes_lookup = gene_lookup_spec.loc[genes]
                genes_lookup_set = set(gene_lookup_spec.loc[genes,'gene_id'].values)  
                
                for seq_record in seq_records:
                    gene_id = seq_record.description 

                    if (gene_id in genes_lookup_set):
                        #find which y1000_id was matched
                        y1000_rlookup = (genes_lookup['gene_id'] == gene_id)
                        for gene, tf in y1000_rlookup.items(): 
                            if tf:
                                y1000_id=gene
                        protein_seq = seq_record.seq
                        f.write('>' + genome_name + '_' + gene_id + ' y1000_id=' + y1000_id + ' gene_id=' + gene_id +'\n')
                        f.write(str(protein_seq) + '\n')  #I wonder why some of the bases were in lower case

            else:   #if not(genome_name in missing_specs):    could also include a missing species list
                #load gene_id map based on the species
                gene_lookup_spec_fname = y1000plus_dir + os.path.normpath("y1000plus_tools_data/y1000plus/id_lookups/" + genome_name + '.csv')
                gene_lookup_spec = pd.read_csv(gene_lookup_spec_fname, index_col='y1000_id')

                #Extract peptide sequences from peptide fasta from genome
                protein_dir = os.path.normpath(y1000plus_dir + 'shen_2018_data/0_332yeast_genomes/332_genome_annotations/pep') + os.sep 

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

def extract_cds_seqs(og_genes, fname, y1000_species_subset): 
    #Looks up coding sequences for given list of orthogroup genes 
    #
    ## Does not work for outgroup species

    #os.mkdir(y1000plus_dir + os.path.normpath('y1000plus_tools_data/y1000plus/proteins_og/' + fname))
    #proteins_og_fname = y1000plus_dir + os.path.normpath('y1000plus_tools_data/y1000plus/proteins_og/' + fname + '/' + fname + '.fasta')

    genome_name_lookup = dict(zip(y1000_species_subset['spec_og_id'],y1000_species_subset['original_genome_id']))

    with open(fname,'w') as f: 
        # group by species
        og_genes_specs = {spec_gene.split('_')[0]: [] for spec_gene in og_genes}
        for y1000_id in og_genes: 
            spec_og_id = y1000_id.split('_')[0]
            og_genes_specs[spec_og_id].append(y1000_id)

        for spec_og_id, genes in og_genes_specs.items(): 
            genome_name = genome_name_lookup[int(spec_og_id)]
            print(genome_name)
            #if S.Cer or C. Albicans, do slightly different routine

            if genome_name == 'saccharomyces_cerevisiae':
                #load gene_id map based on the species
                gene_lookup_spec_fname = y1000plus_dir + os.path.normpath("y1000plus_tools_data/y1000plus/id_lookups/" + genome_name + '.csv')
                gene_lookup_spec = pd.read_csv(gene_lookup_spec_fname, index_col='y1000_id')

                cds_fname = y1000plus_dir + os.path.normpath('shen_2018_data/0_332yeast_genomes/332_genome_annotations/Saccharomyces_cerevisiae_S288C_coding.fasta')

                seq_records = SeqIO.parse(cds_fname, "fasta")

                genes_lookup = gene_lookup_spec.loc[genes]
                genes_lookup_set = set(gene_lookup_spec.loc[genes,'gene_id'].values)  #SC specific

                for seq_record in seq_records:
                    gene_id = seq_record.description.split()[0] #SC specific
                    #print(gene_full)
                    if (gene_id in genes_lookup_set):
                        #find which y1000_id was matched
                        y1000_rlookup = (genes_lookup['gene_id'] == gene_id)
                        for gene, tf in y1000_rlookup.items(): 
                            if tf:
                                y1000_id=gene
                        cds_seq = seq_record.seq
                        f.write('>' + genome_name + '_' + gene_id + ' y1000_id=' + y1000_id + ' gene_id=' + gene_id +
                                ' common_name=' + seq_record.description.split()[1] + #this adds in the cds from the original description
                                '\n')
                        f.write(str(cds_seq) + '\n')  
            elif genome_name == 'candida_albicans':
                #load gene_id map based on the species
                gene_lookup_spec_fname = y1000plus_dir + os.path.normpath("y1000plus_tools_data/y1000plus/id_lookups/" + genome_name + '.csv')
                gene_lookup_spec = pd.read_csv(gene_lookup_spec_fname, index_col='y1000_id')

                cds_fname = y1000plus_dir + os.path.normpath('shen_2018_data/0_332yeast_genomes/332_genome_annotations/Candida_albicans_SC5314_A22_current_default_coding.fasta')
                seq_records = SeqIO.parse(cds_fname, "fasta")

                genes_lookup = gene_lookup_spec.loc[genes]
                genes_lookup_set = set(gene_lookup_spec.loc[genes,'gene_id'].values)  

                for seq_record in seq_records:
                    gene_id = seq_record.description 

                    if (gene_id in genes_lookup_set):
                        #find which y1000_id was matched
                        y1000_rlookup = (genes_lookup['gene_id'] == gene_id)
                        for gene, tf in y1000_rlookup.items(): 
                            if tf:
                                y1000_id=gene
                        cds_seq = seq_record.seq
                        f.write('>' + genome_name + '_' + gene_id + ' y1000_id=' + y1000_id + ' gene_id=' + gene_id +'\n')
                        f.write(str(cds_seq) + '\n')  #I wonder why some of the bases were in lower case

            else:   #if not(genome_name in missing_specs):    could also include a missing species list
                #load gene_id map based on the species
                gene_lookup_spec_fname = y1000plus_dir + os.path.normpath("y1000plus_tools_data/y1000plus/id_lookups/" + genome_name + '.csv')
                gene_lookup_spec = pd.read_csv(gene_lookup_spec_fname, index_col='y1000_id')

                #Extract peptide sequences from peptide fasta from genome
                cds_dir = os.path.normpath(y1000plus_dir + 'shen_2018_data/0_332yeast_genomes/332_genome_annotations/cds') + os.sep 

                cds_fname = cds_dir + genome_name + '.max.cds'

                seq_records = SeqIO.parse(cds_fname, "fasta")

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
                        cds_seq = seq_record.seq
                        f.write('>' + genome_name + '_' + gene_id + ' y1000_id=' + y1000_id + ' gene_full=' + gene_full +
                                ' ' + seq_record.description.split()[2] + #this adds in the cds from the original description
                                '\n')
                        f.write(str(cds_seq) + '\n')  #I wonder why some of the bases were in lower case

    return

def seq_key_func(seq): 
    #returns a gene_id string for a given seq object generated from a promoter fasta file
    #one easy option could be gene_id = seq.id, but the id is hidden in the description. 
    description_dict = {item.split('=')[0] : item.split('=')[1] for item in seq.description.split()}
    gene_id = description_dict['species'] + '@' + description_dict['gene_id'] 
    
    return gene_id


def plot_tree_proms(goi_pair, prom_phyls, t, y1000_species_subset, proms, motif_names, branch_labels):
    #makes tree ready to render for a goi pair and given promoters. 
    #
    #branch_labels can be 
    # 'all':  puts branch length on top, bootstrap/alrt on the bottom or 
    # 'bootstrap': Just puts bootstrap on top
    
    sacc_families = {'Candida': 'Post_WGH',
                 'Kazachstania': 'Post_WGH',
                 'Nakaseomyces': 'Post_WGH',
                 'Naumovozyma': 'Post_WGH',
                 'Saccharomyces': 'Post_WGH',
                 'Tetrapisispora': 'Post_WGH',
                 'Vanderwaltozyma': 'Post_WGH',
                 'Yueomyces': 'Post_WGH',
                 'Zygosaccharomyces': 'ZT',
                 'Zygotorulaspora': 'ZT',
                 'Torulaspora': 'ZT',
                 'Kluyveromyces': 'KLE',
                 'Lachancea': 'KLE',
                 'Eremothecium': 'KLE',
                 'Ashbya': 'KLE'
                }

    #Color Node by species: 
    sacc_colors = {'KLE': "#deb9f6", #e4cee4",#"#C6AFE9", 
                   'ZT': "YellowGreen",
                   'Post_WGH': "LightYellow" #White" # "LightYellow"
                  }

    post_WGH_colors = {'low':  '#8cc3f6', # '#d3d3fe', #'#3192ff',#'#7eeaf7', ##2DD7ED',      #'#e6fcff', 
                       'high': '#fcbba1'} #'#59E3EB'}  #'#ffebe6'}


    #     node_color_dict = {'KLE': "#C6AFE9",
    #                        'ZT': "YellowGreen",
    #                        'Post_WGH': "LightYellow",   #default color for post WGH
    #                        'low': '#7eeaf7',    #syntenic orthologs of low LFC ohnolog
    #                        'high': '#fcbba1',    #syntenic orthologs of high LFC ohnolog
    #                        'outgroup': 'LightGrey'
    #                       }



    #Load Tree
    #t = Tree(fname_tree, format=1)
    ts = TreeStyle()
    ts.show_leaf_name = True #False
    if branch_labels == 'all':
        ts.show_branch_length = True
    else:
        ts.show_branch_length = False

    # #assign Colors, show support values
    # for node in t.traverse():
    #     if node.is_leaf():
    #         #color node by major clade / family if in Sacch clade
    #         row = y1000_species_subset[y1000_species_subset['original_genome_id']==species]
    #         maj_clade = row['Major clade'].values[0]

    #         if maj_clade == 'Saccharomycetaceae':
    #             genus = row['Genus'].values[0]
    #             node_color = sacc_colors[sacc_families[genus]]
    #             if node.name in prom_phyls[goi_pair]['low']: 
    #                 node_color = post_WGH_colors['low']
    #             elif node.name in prom_phyls[goi_pair]['high']: 
    #                 node_color = post_WGH_colors['high']    
    #         #species == outgroup_orig_genome:
    #         #elif species == 'hanseniaspora_vinae':
    #         #    node_color = 'LightGrey'
    #         else:
    #             node_color = 'LightGrey'
    #             #node_color = maj_clade_colors[maj_clade]

    #         nstyle = NodeStyle()
    #         nstyle['bgcolor']=node_color
    #         node.set_style(nstyle)

    #     else:    # If node is not a leaf, add the support label
    #         node_label = TextFace(node.name)
    #         node.add_face(node_label, column=1, position = "branch-bottom")
    #t.render('%%inline', tree_style=ts)

    #Plots promoters for a given tree. 


    # protein_fname_base = y1000plus_dir + os.path.normpath('proteins_og/' + goi_pair + '_' + og)
    # tree_fname = protein_fname_base + '_aln_trimmed.fasta.treefile'
    t.ladderize()
    L_prom = 700
    height = 15
    seq = '-'*L_prom

    motif_colors = {'PDS': 'yellow', 'TATA': 'blue', 'STRE': 'red'}
    motif_lengths = {'PDS': 3*6, 'TATA': 3*8, 'STRE': 3*5 }  #They are triple the size

    #box params:
    width_box = 40
    height_box = 55

    cmap_STRE = cm.get_cmap('Reds')
    vmin = 0.0
    vmax = 8.0
    cmap_STRE_norm = colors.Normalize(vmin=vmin, vmax=vmax)

    cmap_TATA = cm.get_cmap("Blues")

    # To get rid of a set of species for the visualization
    # if less_nonsacc: 
    #     nonsacc_visualization_subset = pd.read_csv(y1000plus_dir + 'species_visualization_subset.csv')
    #     species_subset = ( set(nonsacc_visualization_subset['original_genome_id']) | \
    #                        set(y1000_species[y1000_species['Major clade']=='Saccharomycetaceae']['original_genome_id']) ) # | \
    #                        #set(y1000_species[y1000_species['Species name']==outgroup]['original_genome_id'])
    #                      #)
    #     y1000_species_subset = y1000_species[(y1000_species['original_genome_id'].isin(species_subset))]
    #     node_subset = []
    #     #For each node in the tree:
    #     for node in t.iter_leaves():  
    #         #Get the promoter sequence with motif info, make it into a motif list
    #         if 'saccharomyces_cerevisiae' in node.name:
    #             species='saccharomyces_cerevisiae'
    #         else: 
    #             species = '_'.join(node.name.split('_')[:-2])
    #         if species in species_subset: 
    #             node_subset.append(node.name)

    #     t.prune(node_subset)
    #     t.ladderize()

    # for node in t.iter_leaves():  
    #     print(node.name)

    #For each node in the tree:
    for node in t.traverse(): 
        #name_face = AttrFace("name",fsize=15)
        #node.add_face(name_face, column=0, position="branch-right")     
        if node.is_leaf():#Get the promoter sequence with motif info, make it into a motif list
            if 'saccharomyces_cerevisiae' in node.name:       
                species='saccharomyces_cerevisiae'
                gene_id = species + '@' + node.name.split(species +'_')[1]
            elif 'candida_albicans' in node.name: 
                species = 'candida_albicans'
                gene_id = species + '@' + node.name.split(species +'_')[1]
            else: 
                species = '_'.join(node.name.split('_')[:-2])
                gene_id = species + '@' + '_'.join(node.name.split('_')[-2:])

            #color node by major clade / family if in Sacch clade
            row = y1000_species_subset[y1000_species_subset['original_genome_id']==species]
            maj_clade = row['Major clade'].values[0]

            if maj_clade == 'Saccharomycetaceae':
                genus = row['Genus'].values[0]
                node_color = sacc_colors[sacc_families[genus]]
                if node.name in prom_phyls[goi_pair]['low']: 
                    node_color = post_WGH_colors['low']
                elif node.name in prom_phyls[goi_pair]['high']: 
                    node_color = post_WGH_colors['high']    
            #species == outgroup_orig_genome:
            #elif species == 'hanseniaspora_vinae':
            #    node_color = 'LightGrey'
            else:
                node_color = 'LightGrey'
                #node_color = maj_clade_colors[maj_clade]

            nstyle = NodeStyle()
            nstyle['bgcolor']=node_color
            node.set_style(nstyle)

            prom_results = proms.loc[gene_id]

        #     simple_motifs = [
        #             # seq.start, seq.end, shape, width, height, fgcolor, bgcolor
        #             [10, 60, ">", None, 10, "black", "red", None],
        #             [120, 150, "<", None, 10, "black", "red", None]
        #     ]
            motifs = []

            for motif_name in motif_names: #, 'PDS']:   #Leaving out PDS
                motif_len = motif_lengths[motif_name]
                motif_color = motif_colors[motif_name]
                if prom_results[motif_name + '_count'] >0:

                    for motif in prom_results[motif_name + '_full_features']:
                        loc = motif[0]
                        if loc <= L_prom:
                            direction = motif[1]
                            shape = '>'
                            start = L_prom-loc
                            stop = L_prom-loc + motif_len
                            if direction == 'rev':
                                shape = '<'
                                start = L_prom-loc-motif_len
                                stop = L_prom -loc
                            motifs.append([start,stop,shape,None, height, "black", motif_color, None])

            seqFace = SeqMotifFace(seq, motifs=motifs, seq_format="-")
            node.add_face(seqFace, column=0, position="aligned")


            #Add face for STRE count within 700bp

            motif_name = 'STRE'
            L_STRE_count = 700

            N_STRE = 0
            for result in prom_results['STRE_full_features']:
                if result[0]<L_STRE_count:
                    N_STRE = N_STRE+1

            rgb = colors.to_hex(cmap_STRE(cmap_STRE_norm(N_STRE)))

            rectFace_STRE = RectFace(width=width_box, height=height_box, fgcolor='black', bgcolor=rgb, 
                                label= {"text": str(N_STRE), 
                                        "color": 'black', 
                                        "fontsize": 20, 
                                        "font": 'Verdana'
                                       }
                               )
            node.add_face(rectFace_STRE, column = 1, position= "aligned")

            motif_name = 'TATA'
            L_TATA = 300

            TATA = 0.0
            for result in prom_results['TATA_full_features']:
                if result[0]<L_TATA:
                    TATA = 1.0

            rgb = colors.to_hex(cmap_TATA(TATA))

            rectFace_TATA = RectFace(width=width_box, height=height_box, fgcolor='black', bgcolor=rgb)

            node.add_face(rectFace_TATA, column=2, position="aligned")
        else:  # If node is not a leaf, add the support label
            if branch_labels=='all': 
                node_label = TextFace(node.name)
                node.add_face(node_label, column=1, position = "branch-bottom")
            elif branch_labels =='bootstrap':
                node_label = TextFace(node.name.split('/')[0])
                node.add_face(node_label, column=1, position = "branch-top")
            else: 
                raise ValueError('invalid value for branch_labels: {}'.format(branch_labels))


    return t, ts

def plot_tree_boot_alrt(goi_pair, prom_phyls, fname_tree, fname_out, y1000_species_subset): 

    #og = goi_pair_og_lookup[goi_pair]

    sacc_families = {'Candida': 'Post_WGH',
                     'Kazachstania': 'Post_WGH',
                     'Nakaseomyces': 'Post_WGH',
                     'Naumovozyma': 'Post_WGH',
                     'Saccharomyces': 'Post_WGH',
                     'Tetrapisispora': 'Post_WGH',
                     'Vanderwaltozyma': 'Post_WGH',
                     'Yueomyces': 'Post_WGH',
                     'Zygosaccharomyces': 'ZT',
                     'Zygotorulaspora': 'ZT',
                     'Torulaspora': 'ZT',
                     'Kluyveromyces': 'KLE',
                     'Lachancea': 'KLE',
                     'Eremothecium': 'KLE',
                     'Ashbya': 'KLE'
                    }

    #Color Node by species: 
    sacc_colors = {'KLE': "#deb9f6", #e4cee4",#"#C6AFE9", 
                   'ZT': "YellowGreen",
                   'Post_WGH': "LightYellow" #White" # "LightYellow"
                  }

    post_WGH_colors = {'low':  '#8cc3f6', # '#d3d3fe', #'#3192ff',#'#7eeaf7', ##2DD7ED',      #'#e6fcff', 
                       'high': '#fcbba1'} #'#59E3EB'}  #'#ffebe6'}


    t = Tree(fname_tree, format=1)
    t.ladderize()
    ts = TreeStyle()
    ts.show_leaf_name = True
    ts.show_branch_length = True

    for node in t.traverse():
        if node.is_leaf():
            #Get the species name
            if 'saccharomyces_cerevisiae' in node.name:       
                species='saccharomyces_cerevisiae'
                #gene_id = species + '@' + node.name.split(species +'_')[1]
            elif 'candida_albicans' in node.name: 
                species = 'candida_albicans'
                #gene_id = species + '@' + node.name.split(species +'_')[1]
            else: 
                species = '_'.join(node.name.split('_')[:-2])
                #gene_id = species + '@' + '_'.join(node.name.split('_')[-2:])

            #color node by major clade / family if in Sacch clade
            row = y1000_species_subset[y1000_species_subset['original_genome_id']==species]
            maj_clade = row['Major clade'].values[0]

            if maj_clade == 'Saccharomycetaceae':
                genus = row['Genus'].values[0]
                node_color = sacc_colors[sacc_families[genus]]
                if node.name in prom_phyls[goi_pair]['low']: 
                    node_color = post_WGH_colors['low']
                elif node.name in prom_phyls[goi_pair]['high']: 
                    node_color = post_WGH_colors['high']    
            #species == outgroup_orig_genome:
            #elif species == 'hanseniaspora_vinae':
            #    node_color = 'LightGrey'
            else:
                node_color = 'LightGrey'
                #node_color = maj_clade_colors[maj_clade]

            nstyle = NodeStyle()
            nstyle['bgcolor']=node_color
            node.set_style(nstyle)

        else:    # If node is not a leaf, add the support label
            node_label = TextFace(node.name)
            node.add_face(node_label, column=1, position = "branch-bottom")
    t.render(fname_out)
    #t.render('%%inline', tree_style=ts)
    
def prom_scan_example(goi_pair, promoters_fname, y1000_species_subset):
    #Load Promoter file and find STRE/PDS/TATA box in promoters outputs proms for plotting
    
    output_format = 'full'
    motif_dict = {'STRE': ('CCCCT',700), 'TATA': ('TATA[AT]A[AT][AG]',700), 'PDS': ('AGGGAT',700)}   
    sequence_context = 2
    #L_prom = 1000


    #promoters_fnames = os.listdir(data_processing_dir + os.path.normpath('promoter_phylogeny/promoter_sets') )
    #og_promoter_fnames = { fname.split('_')[0] : fname for fname in promoters_fnames }
    #promoters_fname = data_processing_dir + os.path.normpath('promoter_phylogeny/promoter_sets/') + os.sep + og_promoter_fnames[og]

    proms = yeast_esr_exp.exact_promoter_scan_from_fasta(promoters_fname, motif_dict, output_format = output_format, sequence_context = sequence_context, seq_key_func=seq_key_func)

    proms['spec_genome_name'] = [item.split('@')[0] for item in proms.index]

    proms.index.name = 'seq_id'
    proms.reset_index(inplace=True)
    merge_columns = y1000_species_subset.loc[:,['original_genome_id', 'species_names_fig2']].rename(columns={'original_genome_id':'spec_genome_name'})
    proms = proms.merge(merge_columns, how = 'left', on = 'spec_genome_name')
    proms.set_index('seq_id', inplace=True)

    return proms

def plot_tree_summary(goi_pair, fname_tree, y1000_species_subset, anc_nodes, outgroup_subtract_set, prom_phyls, proms):
    #Plots a summary tree given node pairs to summarize along with the average of the STRE count within a range
    #defined in this function as L_STRE_count
    #returns t_abbrev, an ete3 tree object which is a subset of the original tree. 
    
    sacc_families = {'Candida': 'Post_WGH',
                     'Kazachstania': 'Post_WGH',
                     'Nakaseomyces': 'Post_WGH',
                     'Naumovozyma': 'Post_WGH',
                     'Saccharomyces': 'Post_WGH',
                     'Tetrapisispora': 'Post_WGH',
                     'Vanderwaltozyma': 'Post_WGH',
                     'Yueomyces': 'Post_WGH',
                     'Zygosaccharomyces': 'ZT',
                     'Zygotorulaspora': 'ZT',
                     'Torulaspora': 'ZT',
                     'Kluyveromyces': 'KLE',
                     'Lachancea': 'KLE',
                     'Eremothecium': 'KLE',
                     'Ashbya': 'KLE'
                    }

    #Color Node by species: 
    sacc_colors = {'KLE': "#deb9f6", #e4cee4",#"#C6AFE9", 
                   'ZT': "YellowGreen",
                   'Post_WGH': "LightYellow" #White" # "LightYellow"
                  }

    post_WGH_colors = {'low':  '#8cc3f6', # '#d3d3fe', #'#3192ff',#'#7eeaf7', ##2DD7ED',      #'#e6fcff', 
                       'high': '#fcbba1'} #'#59E3EB'}  #'#ffebe6'}


    node_color_dict = {'KLE':"#deb9f6",
                       'ZT': "YellowGreen",
                       'Post_WGH': "LightYellow",   #default color for post WGH
                       'low': '#8cc3f6',    #syntenic orthologs of low LFC ohnolog
                       'high': '#fcbba1',    #syntenic orthologs of high LFC ohnolog
                       'outgroup': 'LightGrey'
                      }

    L_STRE_count = 700
    cmap_STRE = cm.get_cmap('Reds')
    vmin = 0.0
    vmax = 7.0
    cmap_STRE_norm = colors.Normalize(vmin=vmin, vmax=vmax)
    #box params:
    #width_box = 48
    #height_box = 35
    #box params:
    width_box = 48
    height_box = 35
    height_box_summary = 55

    t = Tree(fname_tree, format=1)
    t.ladderize()
    ts = TreeStyle()
    ts.scale =  60 # 120 pixels per branch length unit
    ts.show_leaf_name = False#True
    # ts.show_branch_length = True
    # ts.mode = "c"
    # ts.arc_start = -90 # 0 degrees = 3 o'clock
    # ts.arc_span = 180
    t_abbrev = t.copy()



    for node in t_abbrev.traverse():
        if node.is_leaf():
            #Get the species and gets the gene_id for finding promoter sequences
            if 'saccharomyces_cerevisiae' in node.name:       
                species='saccharomyces_cerevisiae'
                gene_id = species + '@' + node.name.split(species +'_')[1]
            elif 'candida_albicans' in node.name: 
                species = 'candida_albicans'
                gene_id = species + '@' + node.name.split(species +'_')[1]
            else: 
                species = '_'.join(node.name.split('_')[:-2])
                gene_id = species + '@' + '_'.join(node.name.split('_')[-2:])


            #color node by major clade / family if in Sacch clade
            row = y1000_species_subset[y1000_species_subset['original_genome_id']==species]
            maj_clade = row['Major clade'].values[0]

            if maj_clade == 'Saccharomycetaceae':
                genus = row['Genus'].values[0]
                node_color = sacc_colors[sacc_families[genus]]
                if node.name in prom_phyls[goi_pair]['low']: 
                    node_color = post_WGH_colors['low']
                elif node.name in prom_phyls[goi_pair]['high']: 
                    node_color = post_WGH_colors['high']    
            #species == outgroup_orig_genome:
            #elif species == 'hanseniaspora_vinae':
            #    node_color = 'LightGrey'
            else:
                node_color = 'LightGrey'
                #node_color = maj_clade_colors[maj_clade]

            nstyle = NodeStyle()
            nstyle['bgcolor']=node_color
            node.set_style(nstyle)

            ##Add face for node distance
            #node_dist = TextFace('{:.3f}'.format(node.dist), fsize=8, fstyle='italic')
            #node.add_face(node_dist, column=1, position = "branch-bottom")

            #node label
            node_name_label = AttrFace("name", fsize=10, fgcolor="grey")
            node.add_face(node_name_label, column=1, position = "aligned")

            #add face for number of STREs within L_STRE_count basepairs        
            prom_results = proms.loc[gene_id]

            N_STRE = 0
            for result in prom_results['STRE_full_features']:
                if result[0]<L_STRE_count:
                    N_STRE = N_STRE+1

            rgb = colors.to_hex(cmap_STRE(cmap_STRE_norm(N_STRE)))

            rectFace_STRE = RectFace(width=width_box, height=height_box, fgcolor='black', bgcolor=rgb, 
                                label= {"text": str(N_STRE), 
                                        "color": 'grey', 
                                        "fontsize": 10,
                                        "font": 'Verdana'
                                       }
                               )
            node.add_face(rectFace_STRE, column = 2, position= "aligned")



        else:    # If node is not a leaf, add the support label.  Red if support greater than 90
            
            node_val=node.name.split('/')[0]
            
            if node_val != '':
                if float(node_val)>95:
                    support_color = 'black'
                elif float(node_val)>80:
                    support_color = 'blue'
                else: 
                    support_color = 'red'
                node_label = TextFace(' ' + node.name.split('/')[0], fsize=9, fgcolor = support_color)
                node.add_face(node_label, column=1, position = "branch-bottom")
                #node_dist = TextFace('{:.3f}'.format(node.dist), fsize=8, fstyle='italic')
                #node.add_face(node_dist, column=1, position = "branch-bottom")


    left_leaves = set(t_abbrev.get_leaf_names())
    for node_name, (node1, node2, node_group) in anc_nodes[goi_pair].items(): 
        anc_node = t_abbrev.get_common_ancestor([node1, node2])
        #print(anc_node.name)
        if len(anc_node.name.split('/'))==2:
            anc_node.add_features(shalrt=anc_node.name.split('/')[0], bootstrap=anc_node.name.split('/')[1])
        else:
            anc_node.add_features(shalrt=np.nan,bootstrap = np.nan)
        anc_node.add_features()
        anc_node.name = node_name
        left_leaves = left_leaves - set(anc_node.get_leaf_names())
        nstyle = NodeStyle()
        node_color = node_color_dict[node_group]
        nstyle['bgcolor']=node_color
        anc_node.set_style(nstyle)

        #grouped node label
        summarized_nodes = anc_node.get_leaves()
        
        grouped_node_face = TextFace(anc_node.name + " ({})".format(len(summarized_nodes)), fsize=16, fgcolor="black")
        anc_node.add_face(grouped_node_face, column=1, position = "aligned")
        #print(left_leaves)

        #Add a face for STRE numbers      

        gene_ids = []
        for node in summarized_nodes:
            if 'saccharomyces_cerevisiae' in node.name:       
                species='saccharomyces_cerevisiae'
                gene_id = species + '@' + node.name.split(species +'_')[1]
            elif 'candida_albicans' in node.name: 
                species = 'candida_albicans'
                gene_id = species + '@' + node.name.split(species +'_')[1]
            else: 
                species = '_'.join(node.name.split('_')[:-2])
                gene_id = species + '@' + '_'.join(node.name.split('_')[-2:])
            gene_ids.append(gene_id)

        STRE_counts = []
        for gene_id in gene_ids: 
            prom_results = proms.loc[gene_id]
            N_STRE = 0
            for result in prom_results['STRE_full_features']:
                if result[0]<L_STRE_count:
                    N_STRE = N_STRE+1
            STRE_counts.append(N_STRE)

        rgb = colors.to_hex(cmap_STRE(cmap_STRE_norm(np.mean(STRE_counts))))

        rectFace_STRE = RectFace(width=width_box, height=height_box_summary, fgcolor='black', bgcolor=rgb, 
                                label= {"text": '{:0.2f}'.format(np.mean(STRE_counts)), 
                                        "color": 'black', 
                                        "fontsize": 12, 
                                        "font": 'Verdana'
                                       }
                                 )
        anc_node.add_face(rectFace_STRE, column = 2, position= "aligned")


    left_leaves = left_leaves - outgroup_subtract_set[goi_pair]

    t_abbrev.prune(list(anc_nodes[goi_pair].keys()) + list(left_leaves))

    return t_abbrev, ts


def plot_tree_proms_exp_data(goi_pair, prom_phyls, t, y1000_species_subset, proms, motif_names, branch_labels, goi_exp_data, exp_subset_goi):

    #makes tree ready to render for a goi pair and given promoters. 
    #Includes visualization of gene expression data
    #
    #branch_labels can be 
    # 'all':  puts branch length on top, bootstrap/alrt on the bottom or 
    # 'bootstrap': Just puts bootstrap on top

    #Color Node by species: 
    sacc_colors = {'KLE': "#deb9f6", #e4cee4",#"#C6AFE9", 
                   'ZT': "YellowGreen",
                   'Post_WGH': "LightYellow" #White" # "LightYellow"
                  }

    post_WGH_colors = {'low':  '#8cc3f6', # '#d3d3fe', #'#3192ff',#'#7eeaf7', ##2DD7ED',      #'#e6fcff', 
                       'high': '#fcbba1'} #'#59E3EB'}  #'#ffebe6'}


    ts = TreeStyle()
    ts.show_leaf_name = True #False
    if branch_labels == 'all':
        ts.show_branch_length = True
    else:
        ts.show_branch_length = False

    t.ladderize()
    L_prom = 700
    height = 15
    seq = '-'*L_prom

    motif_colors = {'PDS': 'yellow', 'TATA': 'blue', 'STRE': 'red'}
    motif_lengths = {'PDS': 3*6, 'TATA': 3*8, 'STRE': 3*5 }  #They are triple the size

    #box params:
    width_box = 40
    height_box = 55

    cmap_STRE = cm.get_cmap('Reds')
    vmin = 0.0
    vmax = 8.0
    cmap_STRE_norm = colors.Normalize(vmin=vmin, vmax=vmax)

    cmap_TATA = cm.get_cmap("Blues")

    cmap_exp = cm.get_cmap('coolwarm')
    vmin = -7.5
    vmax = 7.5
    cmap_exp_norm = colors.Normalize(vmin=vmin, vmax=vmax)
    nan_color = '#808080'

    conds = list(goi_exp_data.columns)


    # To get rid of a set of species for the visualization
    # if less_nonsacc: 
    #     nonsacc_visualization_subset = pd.read_csv(y1000plus_dir + 'species_visualization_subset.csv')
    #     species_subset = ( set(nonsacc_visualization_subset['original_genome_id']) | \
    #                        set(y1000_species[y1000_species['Major clade']=='Saccharomycetaceae']['original_genome_id']) ) # | \
    #                        #set(y1000_species[y1000_species['Species name']==outgroup]['original_genome_id'])
    #                      #)
    #     y1000_species_subset = y1000_species[(y1000_species['original_genome_id'].isin(species_subset))]
    #     node_subset = []
    #     #For each node in the tree:
    #     for node in t.iter_leaves():  
    #         #Get the promoter sequence with motif info, make it into a motif list
    #         if 'saccharomyces_cerevisiae' in node.name:
    #             species='saccharomyces_cerevisiae'
    #         else: 
    #             species = '_'.join(node.name.split('_')[:-2])
    #         if species in species_subset: 
    #             node_subset.append(node.name)

    #     t.prune(node_subset)
    #     t.ladderize()

    # for node in t.iter_leaves():  
    #     print(node.name)

    #For each node in the tree:
    for node in t.traverse(): 
        #name_face = AttrFace("name",fsize=15)
        #node.add_face(name_face, column=0, position="branch-right")     
        if node.is_leaf():#Get the promoter sequence with motif info, make it into a motif list
            if 'saccharomyces_cerevisiae' in node.name:       
                species='saccharomyces_cerevisiae'
                gene_id = species + '@' + node.name.split(species +'_')[1]
            elif 'candida_albicans' in node.name: 
                species = 'candida_albicans'
                gene_id = species + '@' + node.name.split(species +'_')[1]
            else: 
                species = '_'.join(node.name.split('_')[:-2])
                gene_id = species + '@' + '_'.join(node.name.split('_')[-2:])


            #Color node based on species group

            spec_group = exp_subset_goi[node.name][1]

            if spec_group in {'high', 'low'}:
                node_color = post_WGH_colors[spec_group]
            else:
                node_color = sacc_colors[spec_group]


            nstyle = NodeStyle()
            nstyle['bgcolor']=node_color
            node.set_style(nstyle)

            prom_results = proms.loc[gene_id]

        #     simple_motifs = [
        #             # seq.start, seq.end, shape, width, height, fgcolor, bgcolor
        #             [10, 60, ">", None, 10, "black", "red", None],
        #             [120, 150, "<", None, 10, "black", "red", None]
        #     ]
            motifs = []

            for motif_name in motif_names: #, 'PDS']:   #Leaving out PDS
                motif_len = motif_lengths[motif_name]
                motif_color = motif_colors[motif_name]
                if prom_results[motif_name + '_count'] >0:

                    for motif in prom_results[motif_name + '_full_features']:
                        loc = motif[0]
                        if loc <= L_prom:
                            direction = motif[1]
                            shape = '>'
                            start = L_prom-loc
                            stop = L_prom-loc + motif_len
                            if direction == 'rev':
                                shape = '<'
                                start = L_prom-loc-motif_len
                                stop = L_prom -loc
                            motifs.append([start,stop,shape,None, height, "black", motif_color, None])

            seqFace = SeqMotifFace(seq, motifs=motifs, seq_format="-")
            node.add_face(seqFace, column=0, position="aligned")


            #Add face for STRE count within 700bp

            motif_name = 'STRE'
            L_STRE_count = 700

            N_STRE = 0
            for result in prom_results['STRE_full_features']:
                if result[0]<L_STRE_count:
                    N_STRE = N_STRE+1

            rgb = colors.to_hex(cmap_STRE(cmap_STRE_norm(N_STRE)))

            rectFace_STRE = RectFace(width=width_box, height=height_box, fgcolor='black', bgcolor=rgb, 
                                label= {"text": str(N_STRE), 
                                        "color": 'black', 
                                        "fontsize": 20, 
                                        "font": 'Verdana'
                                       }
                               )
            node.add_face(rectFace_STRE, column = 1, position= "aligned")

            #Add face for TATA box within 700bp
            motif_name = 'TATA'
            L_TATA = 300

            TATA = 0.0
            for result in prom_results['TATA_full_features']:
                if result[0]<L_TATA:
                    TATA = 1.0

            rgb = colors.to_hex(cmap_TATA(TATA))

            rectFace_TATA = RectFace(width=width_box, height=height_box, fgcolor='black', bgcolor=rgb)

            node.add_face(rectFace_TATA, column=2, position="aligned")


            #add faces for expression data

            for jj, cond in enumerate(conds):
                data_val = goi_exp_data.loc[exp_subset_goi[node.name][0],cond]

                if np.isnan(data_val): 
                    rgb = nan_color
                else: 
                    rgb = colors.to_hex(cmap_exp(cmap_exp_norm(data_val)))

                rectFace_exp_cond = RectFace(width=width_box, height=height_box, fgcolor='black', bgcolor=rgb)
                node.add_face(rectFace_exp_cond, column = 3+jj, position= "aligned")


        else:  # If node is not a leaf, add the support label
            if branch_labels=='all': 
                node_label = TextFace(node.name)
                node.add_face(node_label, column=1, position = "branch-bottom")
            elif branch_labels =='bootstrap':
                node_label = TextFace(node.name.split('/')[0])
                node.add_face(node_label, column=1, position = "branch-top")
            else: 
                raise ValueError('invalid value for branch_labels: {}'.format(branch_labels))


    return t, ts


def background_from_promoters(all_promoters_fname, L_prom):
    background_counts= {'A':0, 'C':0, 'T':0, 'G':0} #, 'R':0, 'Y':0, 'N':0}
    all_promoters = SeqIO.parse(all_promoters_fname, "fasta")

    for prom in all_promoters:
        start = max(len(prom.seq)-L_prom,0)
        seq_trunc = prom.seq[start:len(prom.seq)]
        for base in 'ACTG': #RYN': 
            newcount = seq_trunc.count(base)
            background_counts[base] = background_counts[base] + newcount


    counts_tot = sum(background_counts.values())
    background = {}
    for base, count in background_counts.items(): 
        background[base] = count/counts_tot

    background_rev = {}
    background_rev['A']=background['T']
    background_rev['C']=background['G']
    background_rev['T']=background['A']
    background_rev['G']=background['C']
    
    return background, background_rev


def promoters_all_max_scores(all_promoters_fname, motif_name, thresh, motif_fname, pseudocount, background, background_rev, L_prom, L_prom_min):
    #calculate all max scores for a given set of promoters given a motif and relevant parameters
       
    #Load Motif
    with open(motif_fname) as f: 
        motif = motifs.read(f, "jaspar")

    motif.pseudocounts = pseudocount
    motif.background = background

    motif_rev = copy.deepcopy(motif)
    motif_rev.background = background_rev


    #iterate through promoters and save max score
    all_promoters = SeqIO.parse(all_promoters_fname, "fasta")
    all_max_scores_dict = {'id':[],'len':[],'max_score': []}

    for prom in all_promoters:
        all_max_scores_dict['id'].append(prom.id)

        prom.seq.alphabet = IUPAC.unambiguous_dna

        #truncate by promoter length
        start = max(len(prom.seq)-L_prom,0)
        seq_trunc = prom.seq[start:len(prom.seq)]

        all_max_scores_dict['len'].append(len(seq_trunc))

        if len(seq_trunc)<= L_prom_min: #motif.length:
            all_max_scores_dict['max_score'].append(np.nan)

        else: 
            #calculate forward scores
            scores = motif.pssm.calculate(seq_trunc)
            max_score_fwd = max(scores)

            #calculate reverse scores
            seq_trunc_rev = seq_trunc.reverse_complement()
            scores_rev = motif_rev.pssm.calculate(seq_trunc_rev)
            max_score_rev = max(scores_rev)

            max_score = max(max_score_fwd,max_score_rev)

            all_max_scores_dict['max_score'].append(max_score)

            ##catch when log odds is above a certain threshold

    all_max_scores_motif = pd.DataFrame.from_dict(all_max_scores_dict, orient='columns')
    all_max_scores_motif.set_index('id', inplace=True)
    ecdf = ECDF(all_max_scores_motif['max_score'])
    all_max_scores_motif['percentile'] = ecdf(all_max_scores_motif['max_score'])
    all_max_scores_motif.rename(columns = {col : motif_name + '_' + col for col in ['max_score', 'percentile']}, inplace = True)
    
    return all_max_scores_motif


def find_spec_synteny(synteny_window_new_cols, focus_spec, focus_anc, genome_name_lookup, y1000_species_subset, anc_genes_ogs_flat, og_genes_lookup, bp_range):
    #Given a species, find genes surrounding the ortholog of the goi pair.


    #Get focus species genes that are present in focus ogs
    spec_og_id_lookup = dict(zip( y1000_species_subset['original_genome_id'], y1000_species_subset['spec_og_id']))
    spec_og_id = spec_og_id_lookup[focus_spec]
    focus_ogs = anc_genes_ogs_flat[focus_anc]

    og_genes = set()
    for og in focus_ogs: 
        og_genes = og_genes | og_genes_lookup[og]

    focus_spec_genes = []
    for og_gene in list(og_genes): 
        spec_og_id_test = int(og_gene.split('_')[0])
        if spec_og_id == spec_og_id_test:
            focus_spec_genes.append(og_gene)

    #Make full_gene lookup for focus spec
    focus_spec_lookup_fname = y1000plus_dir + os.path.normpath('id_lookups/' + focus_spec + '.csv')
    focus_spec_lookup = pd.read_csv(focus_spec_lookup_fname, index_col=0)
    y1000_id_gene_full_lookup = dict(zip(focus_spec_lookup['y1000_id'],focus_spec_lookup['gene_full']))        
    gene_full_y1000_id_lookup = dict(zip(focus_spec_lookup['gene_full'],focus_spec_lookup['y1000_id']))    
    y1000_id_gene_id_lookup = dict(zip(focus_spec_lookup['y1000_id'],focus_spec_lookup.index))    

    # Load GTF for given sequence: 
    gtf_dir = y1000plus_dir + "0_332yeast_genomes/332_genome_annotations/gtf/"
    db_fname = gtf_dir + 'gffutils_dbs/' + focus_spec + '.db'

    gtf_db = gffutils.FeatureDB(db_fname)


    # Make a dictionary of scaffold sizes
    genome_assembly_fname = y1000plus_dir + os.path.normpath('0_332yeast_genomes/332_genome_assemblies/' + focus_spec + '.fas')

    scaffold_sizes = {}
    for seq in SeqIO.parse(genome_assembly_fname, "fasta"):
        scaffold_sizes[seq.id] = len(seq.seq)

    for og_gene in focus_spec_genes: 
        #print(og_gene)
        focus_gene_column = y1000_find_surrounding_genes(og_gene, y1000_id_gene_id_lookup, y1000_id_gene_full_lookup, gene_full_y1000_id_lookup, og_genes_lookup, gtf_db, bp_range, scaffold_sizes, anc_genes_ogs_flat, focus_spec_genes, focus_spec)
        synteny_window_new_cols = synteny_window_new_cols.merge(focus_gene_column, left_on='Ancestor', right_index = True, how='left')

    return synteny_window_new_cols


def y1000_find_surrounding_genes(og_gene, y1000_id_gene_id_lookup, y1000_id_gene_full_lookup, gene_full_y1000_id_lookup, og_genes_lookup, gtf_db, bp_range, scaffold_sizes, anc_genes_ogs_flat, focus_spec_genes, focus_spec):
    #for a given gene find the surrounding genes 
    
    focus_gene_id = y1000_id_gene_id_lookup[og_gene]

    gene_full = y1000_id_gene_full_lookup[og_gene]

    cursor = gtf_db.execute('select * from features where attributes like "%' + gene_full + '%"')
    all_features = cursor.fetchall()
    if len(all_features) == 0:
        print('No features found ' + gene_full + ' ' + focus_spec)  

    CDS_list = [feature for feature in all_features if feature['featuretype']=='CDS']
    if len(CDS_list)==1: 
        start = CDS_list[0]['start']
        end = CDS_list[0]['end']
    elif len(CDS_list)>1: 
        print('More than one CDS - introns?.  N of CDS: '+ str(len(CDS_list)))
        start = min([CDS['start'] for CDS in CDS_list])
        end = max([CDS['end'] for CDS in CDS_list])
    else: 
        raise AssertionError('No CDS found' + og_gene)

    CDS = CDS_list[0]
    strand = CDS['strand']
    scaffold = CDS['seqid']

    # if strand=='+': 
    #     lower = start
    #     higher = end
    # elif strand == '-': 
    #     lower = end
    #     higher = start

    low = max(0,start-bp_range/2)
    high = min(scaffold_sizes[scaffold], end + bp_range/2)

    surrounding_CDS = list(gtf_db.region(region=(scaffold,low,high), completely_within=False, featuretype='CDS'))

    #surrounding_CDS info: (gene_full, gene_id, strand, start, end)
    surrounding_CDS_data = []
    for CDS_surr in surrounding_CDS:
        surrounding_CDS_data.append((CDS_surr.attributes['gene_id'][0],CDS_surr.id,CDS_surr.strand,CDS_surr.start,CDS_surr.end ))

    #Find OG and line that each gene is in
    anc_focus_spec_genes = OrderedDict([(anc,[]) for anc in anc_genes_ogs_flat.keys()])
    for (gene_full, gene_id, strand, start, end) in surrounding_CDS_data:
        y1000_id = gene_full_y1000_id_lookup[gene_full]
        for anc,og_list in anc_genes_ogs_flat.items():
            for og in og_list: 
                og_genes = og_genes_lookup[og]
                if y1000_id in og_genes:  #could put in a test to make sure this only happens once
                    anc_focus_spec_genes[anc].append(y1000_id_gene_id_lookup[y1000_id])
                    #could make triple dashes if len(anc_focus_spec_genes[anc]) == 0:
            #remove redundancies created by genes with introns
            anc_focus_spec_genes[anc]=list(set(anc_focus_spec_genes[anc]))
    #Merge into synteny window dataframe

    focus_gene_column = pd.DataFrame(pd.Series(anc_focus_spec_genes), columns = [focus_spec + '_' + focus_gene_id])



    return focus_gene_column
